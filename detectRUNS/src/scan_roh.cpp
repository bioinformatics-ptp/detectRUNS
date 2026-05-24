#include "scan_roh.h"
#include <stdexcept>
#include <cstring>
#include <cstdio>
#include <R_ext/Print.h>
extern "C" { void R_FlushConsole(void); }

#ifdef _OPENMP
#include <omp.h>
#endif

// Print a 40-char progress bar to stderr via Rprintf; call only from one thread.
static inline void print_progress(int j, int M) {
    const int pct    = (j + 1) * 100 / M;
    const int filled = pct * 40 / 100;
    char bar[41];
    for (int b = 0; b < 40; ++b) bar[b] = (b < filled) ? '#' : '.';
    bar[40] = '\0';
    Rprintf("\r  Progress: [%s] %3d%%  (%d/%d SNPs)", bar, pct, j + 1, M);
    R_FlushConsole();
}


// ===========================================================================
// CONSECUTIVE METHOD
// ===========================================================================

// Per-individual state for the consecutive state machine.
struct ConsecState {
    int32_t  n_snp;           // SNPs in current run (target + accepted opp/miss)
    int32_t  n_opposite;      // accepted opposite genotypes in run
    int32_t  n_missing;       // accepted missing genotypes in run
    int32_t  run_start_bp;
    int32_t  run_end_bp;
    int32_t  run_start_snp;   // BIM index of first SNP in run (for snp_freq)
    int32_t  run_end_snp;     // BIM index of last  SNP in run
    int32_t  last_bp;         // bp position of last processed SNP (gap calc)
    uint8_t  run_chrom_idx;   // chromosome of current run
    uint8_t  last_chrom_idx;  // 0xFF = sentinel (no previous chromosome)
    uint8_t  in_run;
};

static ConsecState make_consec_state() {
    ConsecState s;
    memset(&s, 0, sizeof(s));
    s.last_chrom_idx = 0xFF;
    return s;
}

// Emit the current run if it meets min_snps and min_length_bp.
static void emit_consec(
    const ConsecState& s, int sample_idx,
    const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    if (s.n_snp < p.min_snps) return;
    const int32_t len = s.run_end_bp - s.run_start_bp;
    if (len < p.min_length_bp) return;

    RohRecord rec;
    rec.sample_idx = sample_idx;
    rec.chrom_idx  = s.run_chrom_idx;
    rec.start_bp   = s.run_start_bp;
    rec.end_bp     = s.run_end_bp;
    rec.n_snps     = s.n_snp;
    rec.n_het      = static_cast<int16_t>(s.n_opposite);
    rec.n_missing  = static_cast<int16_t>(s.n_missing);
    out.push_back(rec);

    sum.n_roh++;
    sum.total_length_bp += len;
    sum.n_snps_in_roh   += s.n_snp;

    for (int j = s.run_start_snp; j <= s.run_end_snp; ++j)
        snp_freq[j]++;
}

// Update state for one SNP (called in BIM order, for one individual).
static void update_consec(
    ConsecState& s, int8_t geno,
    const SnpInfo& snp, int snp_idx,
    int sample_idx, const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    // --- Chromosome boundary ---
    // Reset last_bp to current position so the gap calc below yields 0.
    if (snp.chrom_idx != s.last_chrom_idx) {
        if (s.in_run) {
            emit_consec(s, sample_idx, p, out, sum, snp_freq);
            s.in_run = 0;
        }
        s.last_chrom_idx = snp.chrom_idx;
        s.last_bp        = snp.bp_pos;
    }

    // --- Gap check ---
    if (snp.bp_pos - s.last_bp >= p.max_gap) {
        if (s.in_run) {
            emit_consec(s, sample_idx, p, out, sum, snp_freq);
            s.in_run = 0;
        }
    }

    // --- Genotype ---
    const int8_t opp = static_cast<int8_t>(1 - p.target);

    if (geno == p.target) {
        if (!s.in_run) {
            s.in_run        = 1;
            s.n_snp         = 0;
            s.n_opposite    = 0;
            s.n_missing     = 0;
            s.run_start_bp  = snp.bp_pos;
            s.run_start_snp = snp_idx;
            s.run_chrom_idx = snp.chrom_idx;
        }
        s.n_snp++;
        s.run_end_bp  = snp.bp_pos;
        s.run_end_snp = snp_idx;

    } else if (geno == opp) {
        if (s.in_run) {
            if (s.n_opposite < p.max_opposite) {
                s.n_opposite++;
                s.n_snp++;
                s.run_end_bp  = snp.bp_pos;
                s.run_end_snp = snp_idx;
            } else {
                emit_consec(s, sample_idx, p, out, sum, snp_freq);
                s.in_run = 0;
            }
        }

    } else {
        // GENO_MISSING
        if (s.in_run) {
            if (s.n_missing < p.max_missing) {
                s.n_missing++;
                s.n_snp++;
                s.run_end_bp  = snp.bp_pos;
                s.run_end_snp = snp_idx;
            } else {
                emit_consec(s, sample_idx, p, out, sum, snp_freq);
                s.in_run = 0;
            }
        }
    }

    s.last_bp = snp.bp_pos;
}

// Single pass over all SNP rows; update every individual's state per SNP.
// With OpenMP: persistent parallel region; one thread decodes each BED row
// (omp single), all threads share individual processing (omp for).
// Thread-local record and snp_freq buffers are merged into the output at end.
static void scan_consecutive(
    const BedFile& bed, const BimData& bim, const FamData& fam,
    const ScanParams& params,
    std::vector<RohRecord>& records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>& snp_freq)
{
    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());

#ifdef _OPENMP
    const int actual = (params.n_threads > 0) ? params.n_threads : 1;
#else
    const int actual = 1;
#endif

    std::vector<std::vector<RohRecord>> tl_rec(actual);
    std::vector<std::vector<int32_t>>   tl_freq(actual, std::vector<int32_t>(M, 0));

    std::vector<ConsecState> states(N);
    for (int i = 0; i < N; ++i) states[i] = make_consec_state();

    std::vector<int8_t> row_buf(N);

#ifdef _OPENMP
    #pragma omp parallel num_threads(actual)
    {
        const int tid = omp_get_thread_num();
        for (int j = 0; j < M; ++j) {
            #pragma omp single
            {
                decode_snp_row(bed, j, row_buf.data());
                if (params.verbose && (j * 20 / M != (j - 1) * 20 / M || j == M - 1))
                    print_progress(j, M);
            }

            const SnpInfo& snp = bim.snps[j];
            #pragma omp for schedule(static)
            for (int i = 0; i < N; ++i)
                update_consec(states[i], row_buf[i], snp, j, i,
                              params, tl_rec[tid], summaries[i], tl_freq[tid]);
        }
        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i)
            if (states[i].in_run)
                emit_consec(states[i], i, params, tl_rec[tid], summaries[i], tl_freq[tid]);
    }
    if (params.verbose) Rprintf("\n");
#else
    for (int j = 0; j < M; ++j) {
        decode_snp_row(bed, j, row_buf.data());
        if (params.verbose && (j * 20 / M != (j - 1) * 20 / M || j == M - 1))
            print_progress(j, M);
        const SnpInfo& snp = bim.snps[j];
        for (int i = 0; i < N; ++i)
            update_consec(states[i], row_buf[i], snp, j, i,
                          params, tl_rec[0], summaries[i], tl_freq[0]);
    }
    if (params.verbose) Rprintf("\n");
    for (int i = 0; i < N; ++i)
        if (states[i].in_run)
            emit_consec(states[i], i, params, tl_rec[0], summaries[i], tl_freq[0]);
#endif

    for (int t = 0; t < actual; ++t) {
        records.insert(records.end(), tl_rec[t].begin(), tl_rec[t].end());
        for (int j = 0; j < M; ++j) snp_freq[j] += tl_freq[t][j];
    }
}


// ===========================================================================
// SLIDING WINDOW METHOD — streaming ring-buffer implementation (Phase 5)
//
// Replaces the previous per-chromosome buffer approach.
// Memory: O(N × W) instead of O(N × max_chrom_snps).
//
// Algorithm (Bjelland 2013) streamed with W-SNP decision delay:
//   When SNP j is processed, window w = j-W+1 just became complete.
//   SNP s = w is then fully covered by windows [max(0,s-W+1), s], all
//   already evaluated → decide s immediately.
//   At chromosome end, flush the trailing W-1 undecided SNPs.
// ===========================================================================

// Emit a qualifying run and update summary + snp_freq.
static void emit_sw(
    const IndivStateSW& st,
    int sample_idx, uint8_t chrom_idx,
    const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    if (st.n_snp_run < p.min_snps) return;
    const int32_t len = st.run_end_bp - st.run_start_bp;
    if (len < p.min_length_bp) return;

    RohRecord rec;
    rec.sample_idx = sample_idx;
    rec.chrom_idx  = chrom_idx;
    rec.start_bp   = st.run_start_bp;
    rec.end_bp     = st.run_end_bp;
    rec.n_snps     = st.n_snp_run;
    rec.n_het      = static_cast<int16_t>(st.n_opp_run);
    rec.n_missing  = static_cast<int16_t>(st.n_miss_run);
    out.push_back(rec);

    sum.n_roh++;
    sum.total_length_bp += len;
    sum.n_snps_in_roh   += st.n_snp_run;

    for (int k = st.run_start_snp; k < st.run_start_snp + st.n_snp_run; ++k)
        snp_freq[k]++;
}

// Decide whether within-chrom SNP index s is in ROH.
// [from_win, to_win]: range of window indices covering SNP s (already evaluated).
static void decide_snp_sw(
    IndivStateSW& st,
    int s, int from_win, int to_win,
    int sample_idx, int first_snp_idx, uint8_t chrom_idx,
    const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    int h = 0;
    for (int w = from_win; w <= to_win; ++w)
        h += st.win_pass_buf[w % (MAX_WINDOW * 2)];

    const int   n_cov  = to_win - from_win + 1;
    // >= matches PLINK --homozyg behaviour: a SNP where h/n_cov equals the
    // threshold exactly is included in the run (vs old strict-greater-than).
    // Both sides cast to float for consistency with the PED-path snpInRunCpp.
    const bool  in_roh = (static_cast<float>(h) / n_cov >= static_cast<float>(p.threshold));
    const int8_t opp   = static_cast<int8_t>(1 - p.target);

    if (in_roh) {
        const int32_t cur_bp = st.bp_buf[s % MAX_WINDOW];
        // Run-level gap check (PLINK-compatible): if the distance from the
        // last in-ROH SNP exceeds maxGap, close the current run and restart.
        if (st.in_run && cur_bp - st.run_end_bp > p.max_gap) {
            emit_sw(st, sample_idx, chrom_idx, p, out, sum, snp_freq);
            st.in_run = 0;
        }
        if (!st.in_run) {
            st.in_run        = 1;
            st.run_start_bp  = cur_bp;
            st.run_start_snp = first_snp_idx + s;
            st.n_opp_run     = 0;
            st.n_miss_run    = 0;
            st.n_snp_run     = 0;
        }
        const int8_t g = st.geno_buf[s % MAX_WINDOW];
        if      (g == GENO_MISSING) ++st.n_miss_run;
        else if (g == opp)          ++st.n_opp_run;
        ++st.n_snp_run;
        st.run_end_bp = cur_bp;
    } else {
        if (st.in_run) {
            emit_sw(st, sample_idx, chrom_idx, p, out, sum, snp_freq);
            st.in_run = 0;
        }
    }
}

// Process one SNP for one individual (streaming, within-chromosome).
// Evaluates window w = pos-W+1 when it becomes complete, then decides SNP s = w.
static void update_sw(
    IndivStateSW& st,
    int8_t geno, int32_t bp,
    int sample_idx, int first_snp_idx, uint8_t chrom_idx,
    const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    const int    pos = st.chrom_pos;
    const int    W   = p.window_size;
    const int8_t opp = static_cast<int8_t>(1 - p.target);

    st.geno_buf[pos % MAX_WINDOW] = geno;
    st.bp_buf  [pos % MAX_WINDOW] = bp;

    if (pos >= W - 1) {
        const int w = pos - W + 1;

        // Evaluate window w = [w .. pos] — gap check is now run-level (see decide_snp_sw).
        int8_t win_pass = 0;
        {
            int n_opp = 0, n_miss = 0;
            for (int k = w; k <= pos; ++k) {
                const int8_t g = st.geno_buf[k % MAX_WINDOW];
                if      (g == GENO_MISSING) ++n_miss;
                else if (g == opp)          ++n_opp;
            }
            win_pass = (n_opp <= p.max_opposite && n_miss <= p.max_missing) ? 1 : 0;
        }
        st.win_pass_buf[w % (MAX_WINDOW * 2)] = win_pass;

        // SNP s = w is now fully covered; decide its in-ROH status.
        const int s      = w;
        const int from_w = (s - W + 1 > 0) ? (s - W + 1) : 0;
        decide_snp_sw(st, s, from_w, s,
                      sample_idx, first_snp_idx, chrom_idx,
                      p, out, sum, snp_freq);
    }

    st.chrom_pos++;
}

// Called at chromosome end.  Decides the trailing W-1 SNPs (their last
// covering window was the final window of the chromosome, now complete),
// then closes any open run and resets chrom_pos for the next chromosome.
static void flush_chrom_sw(
    IndivStateSW& st,
    int sample_idx, int first_snp_idx, uint8_t chrom_idx,
    const ScanParams& p,
    std::vector<RohRecord>& out,
    IndivSummary& sum,
    std::vector<int32_t>& snp_freq)
{
    const int C = st.chrom_pos;
    const int W = p.window_size;

    if (C >= W) {
        // Inner loop decided SNPs 0 .. last_win where last_win = C-W.
        // Flush SNPs last_win+1 .. C-1 (W-1 SNPs), all using last_win as to_win.
        const int last_win = C - W;
        for (int s = last_win + 1; s < C; ++s) {
            const int from_w = (s - W + 1 > 0) ? (s - W + 1) : 0;
            decide_snp_sw(st, s, from_w, last_win,
                          sample_idx, first_snp_idx, chrom_idx,
                          p, out, sum, snp_freq);
        }
    }
    // C < W: no complete windows → no ROH possible; nothing to decide.

    if (st.in_run) {
        emit_sw(st, sample_idx, chrom_idx, p, out, sum, snp_freq);
        st.in_run = 0;
    }

    st.chrom_pos = 0;
}

// Single-pass streaming sliding window scan.
// O(N × W) working memory (IndivStateSW ring buffers), one BED pass.
// OpenMP: same persistent-region pattern as scan_consecutive.
// Chromosome-boundary flush and shared-variable update both use omp single
// so all threads see consistent first_snp_idx / cur_chrom values.
static void scan_sliding(
    const BedFile& bed, const BimData& bim, const FamData& fam,
    const ScanParams& params,
    std::vector<RohRecord>& records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>& snp_freq)
{
    if (params.window_size > MAX_WINDOW)
        throw std::runtime_error(
            "window_size exceeds MAX_WINDOW (256); reduce windowSize");

    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());
    if (M == 0) return;

#ifdef _OPENMP
    const int actual = (params.n_threads > 0) ? params.n_threads : 1;
#else
    const int actual = 1;
#endif

    std::vector<std::vector<RohRecord>> tl_rec(actual);
    std::vector<std::vector<int32_t>>   tl_freq(actual, std::vector<int32_t>(M, 0));

    std::vector<IndivStateSW> states(N);
    for (int i = 0; i < N; ++i)
        memset(&states[i], 0, sizeof(IndivStateSW));

    std::vector<int8_t> row_buf(N);
    int     first_snp_idx = 0;
    uint8_t cur_chrom     = bim.snps[0].chrom_idx;

#ifdef _OPENMP
    #pragma omp parallel num_threads(actual)
    {
        const int tid = omp_get_thread_num();
        for (int j = 0; j < M; ++j) {
            #pragma omp single
            {
                decode_snp_row(bed, j, row_buf.data());
                if (params.verbose && (j * 20 / M != (j - 1) * 20 / M || j == M - 1))
                    print_progress(j, M);
            }

            const SnpInfo& snp = bim.snps[j];

            if (snp.chrom_idx != cur_chrom) {
                #pragma omp for schedule(static)
                for (int i = 0; i < N; ++i)
                    flush_chrom_sw(states[i], i, first_snp_idx, cur_chrom,
                                   params, tl_rec[tid], summaries[i], tl_freq[tid]);

                #pragma omp single
                {
                    first_snp_idx = j;
                    cur_chrom     = snp.chrom_idx;
                }
            }

            #pragma omp for schedule(static)
            for (int i = 0; i < N; ++i)
                update_sw(states[i], row_buf[i], snp.bp_pos,
                          i, first_snp_idx, cur_chrom, params,
                          tl_rec[tid], summaries[i], tl_freq[tid]);
        }
        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i)
            flush_chrom_sw(states[i], i, first_snp_idx, cur_chrom,
                           params, tl_rec[tid], summaries[i], tl_freq[tid]);
    }
    if (params.verbose) Rprintf("\n");
#else
    for (int j = 0; j < M; ++j) {
        decode_snp_row(bed, j, row_buf.data());
        if (params.verbose && (j * 20 / M != (j - 1) * 20 / M || j == M - 1))
            print_progress(j, M);
        const SnpInfo& snp = bim.snps[j];

        if (snp.chrom_idx != cur_chrom) {
            for (int i = 0; i < N; ++i)
                flush_chrom_sw(states[i], i, first_snp_idx, cur_chrom,
                               params, tl_rec[0], summaries[i], tl_freq[0]);
            first_snp_idx = j;
            cur_chrom     = snp.chrom_idx;
        }

        for (int i = 0; i < N; ++i)
            update_sw(states[i], row_buf[i], snp.bp_pos,
                      i, first_snp_idx, cur_chrom, params,
                      tl_rec[0], summaries[i], tl_freq[0]);
    }
    if (params.verbose) Rprintf("\n");
    for (int i = 0; i < N; ++i)
        flush_chrom_sw(states[i], i, first_snp_idx, cur_chrom,
                       params, tl_rec[0], summaries[i], tl_freq[0]);
#endif

    for (int t = 0; t < actual; ++t) {
        records.insert(records.end(), tl_rec[t].begin(), tl_rec[t].end());
        for (int j = 0; j < M; ++j) snp_freq[j] += tl_freq[t][j];
    }
}


// ===========================================================================
// INDIVIDUAL-MAJOR BED — Phase 7
//
// When the BED file is stored in individual-major layout (one row per sample,
// M genotypes per row), each individual is completely independent.  This makes
// the parallelism trivial: no shared row_buf, no omp single decode, no
// chromosome-boundary barriers across threads.
//
// Outer loop: individuals   → distributed across threads (omp parallel for)
// Inner loop: SNPs in order → sequential per individual (state machine)
//
// Progress bar is printed in the sequential (single-thread) path only.
// ===========================================================================

#ifndef _OPENMP
static inline void print_progress_indiv(int i, int N) {
    const int pct    = (i + 1) * 100 / N;
    const int filled = pct * 40 / 100;
    char bar[41];
    for (int b = 0; b < 40; ++b) bar[b] = (b < filled) ? '#' : '.';
    bar[40] = '\0';
    Rprintf("\r  Progress: [%s] %3d%%  (%d/%d individuals)", bar, pct, i + 1, N);
    R_FlushConsole();
}
#endif

static void scan_consecutive_indiv_major(
    const BedFile& bed, const BimData& bim, const FamData& fam,
    const ScanParams& params,
    std::vector<RohRecord>& records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>& snp_freq)
{
    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());

#ifdef _OPENMP
    const int actual = (params.n_threads > 0) ? params.n_threads : 1;
#else
    const int actual = 1;
#endif

    std::vector<std::vector<RohRecord>> tl_rec(actual);
    std::vector<std::vector<int32_t>>   tl_freq(actual, std::vector<int32_t>(M, 0));

#ifdef _OPENMP
    #pragma omp parallel num_threads(actual)
    {
        const int tid = omp_get_thread_num();
        std::vector<int8_t> ind_buf(M);

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            decode_individual_row(bed, i, ind_buf.data());
            ConsecState s = make_consec_state();
            for (int j = 0; j < M; ++j)
                update_consec(s, ind_buf[j], bim.snps[j], j, i, params,
                              tl_rec[tid], summaries[i], tl_freq[tid]);
            if (s.in_run)
                emit_consec(s, i, params, tl_rec[tid], summaries[i], tl_freq[tid]);
        }
    }
    if (params.verbose) Rprintf("\n");
#else
    std::vector<int8_t> ind_buf(M);
    for (int i = 0; i < N; ++i) {
        decode_individual_row(bed, i, ind_buf.data());
        if (params.verbose && (i * 20 / N != (i - 1) * 20 / N || i == N - 1))
            print_progress_indiv(i, N);
        ConsecState s = make_consec_state();
        for (int j = 0; j < M; ++j)
            update_consec(s, ind_buf[j], bim.snps[j], j, i, params,
                          tl_rec[0], summaries[i], tl_freq[0]);
        if (s.in_run)
            emit_consec(s, i, params, tl_rec[0], summaries[i], tl_freq[0]);
    }
    if (params.verbose) Rprintf("\n");
#endif

    for (int t = 0; t < actual; ++t) {
        records.insert(records.end(), tl_rec[t].begin(), tl_rec[t].end());
        for (int j = 0; j < M; ++j) snp_freq[j] += tl_freq[t][j];
    }
}

static void scan_sliding_indiv_major(
    const BedFile& bed, const BimData& bim, const FamData& fam,
    const ScanParams& params,
    std::vector<RohRecord>& records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>& snp_freq)
{
    if (params.window_size > MAX_WINDOW)
        throw std::runtime_error(
            "window_size exceeds MAX_WINDOW (256); reduce windowSize");

    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());
    if (M == 0) return;

#ifdef _OPENMP
    const int actual = (params.n_threads > 0) ? params.n_threads : 1;
#else
    const int actual = 1;
#endif

    std::vector<std::vector<RohRecord>> tl_rec(actual);
    std::vector<std::vector<int32_t>>   tl_freq(actual, std::vector<int32_t>(M, 0));

    const uint8_t first_chrom = bim.snps[0].chrom_idx;

#ifdef _OPENMP
    #pragma omp parallel num_threads(actual)
    {
        const int tid = omp_get_thread_num();
        std::vector<int8_t> ind_buf(M);

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            decode_individual_row(bed, i, ind_buf.data());

            IndivStateSW st;
            memset(&st, 0, sizeof(IndivStateSW));
            int     first_snp_idx = 0;
            uint8_t cur_chrom     = first_chrom;

            for (int j = 0; j < M; ++j) {
                const SnpInfo& snp = bim.snps[j];
                if (snp.chrom_idx != cur_chrom) {
                    flush_chrom_sw(st, i, first_snp_idx, cur_chrom, params,
                                   tl_rec[tid], summaries[i], tl_freq[tid]);
                    first_snp_idx = j;
                    cur_chrom     = snp.chrom_idx;
                }
                update_sw(st, ind_buf[j], snp.bp_pos, i, first_snp_idx, cur_chrom,
                          params, tl_rec[tid], summaries[i], tl_freq[tid]);
            }
            flush_chrom_sw(st, i, first_snp_idx, cur_chrom, params,
                           tl_rec[tid], summaries[i], tl_freq[tid]);
        }
    }
    if (params.verbose) Rprintf("\n");
#else
    std::vector<int8_t> ind_buf(M);
    for (int i = 0; i < N; ++i) {
        decode_individual_row(bed, i, ind_buf.data());
        if (params.verbose && (i * 20 / N != (i - 1) * 20 / N || i == N - 1))
            print_progress_indiv(i, N);

        IndivStateSW st;
        memset(&st, 0, sizeof(IndivStateSW));
        int     first_snp_idx = 0;
        uint8_t cur_chrom     = first_chrom;

        for (int j = 0; j < M; ++j) {
            const SnpInfo& snp = bim.snps[j];
            if (snp.chrom_idx != cur_chrom) {
                flush_chrom_sw(st, i, first_snp_idx, cur_chrom, params,
                               tl_rec[0], summaries[i], tl_freq[0]);
                first_snp_idx = j;
                cur_chrom     = snp.chrom_idx;
            }
            update_sw(st, ind_buf[j], snp.bp_pos, i, first_snp_idx, cur_chrom,
                      params, tl_rec[0], summaries[i], tl_freq[0]);
        }
        flush_chrom_sw(st, i, first_snp_idx, cur_chrom, params,
                       tl_rec[0], summaries[i], tl_freq[0]);
    }
    if (params.verbose) Rprintf("\n");
#endif

    for (int t = 0; t < actual; ++t) {
        records.insert(records.end(), tl_rec[t].begin(), tl_rec[t].end());
        for (int j = 0; j < M; ++j) snp_freq[j] += tl_freq[t][j];
    }
}


// ===========================================================================
// PUBLIC ENTRY POINT
// ===========================================================================

void scan_roh_bed(
    const BedFile& bed, const BimData& bim, const FamData& fam,
    const ScanParams& params,
    std::vector<RohRecord>& records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>& snp_freq)
{
    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());

    summaries.assign(N, {0, 0LL, 0});
    snp_freq.assign(M, 0);
    records.clear();

    if (N == 0 || M == 0) return;

    if (bed.layout == BedLayout::SNP_MAJOR) {
        if (params.method == 0)
            scan_consecutive(bed, bim, fam, params, records, summaries, snp_freq);
        else
            scan_sliding(bed, bim, fam, params, records, summaries, snp_freq);
    } else {
        if (params.method == 0)
            scan_consecutive_indiv_major(bed, bim, fam, params, records, summaries, snp_freq);
        else
            scan_sliding_indiv_major(bed, bim, fam, params, records, summaries, snp_freq);
    }
}
