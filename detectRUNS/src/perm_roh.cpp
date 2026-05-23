#include "perm_roh.h"
#include "bed_reader.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <cstring>
#include <stdexcept>
#include <R_ext/Print.h>
extern "C" { void R_FlushConsole(void); }

#ifdef _OPENMP
#include <omp.h>
#endif


// ===========================================================================
// Helpers
// ===========================================================================

// Build the list of per-chromosome SNP ranges in BIM order.
static std::vector<ChromRange> build_chrom_ranges(const BimData& bim)
{
    std::vector<ChromRange> ranges;
    const int M = static_cast<int>(bim.snps.size());
    if (M == 0) return ranges;

    ChromRange cur;
    cur.chrom_idx = bim.snps[0].chrom_idx;
    cur.first_snp = 0;

    for (int j = 1; j <= M; ++j) {
        bool boundary = (j == M) || (bim.snps[j].chrom_idx != cur.chrom_idx);
        if (boundary) {
            cur.n_snps = j - cur.first_snp;
            ranges.push_back(cur);
            if (j < M) {
                cur.chrom_idx = bim.snps[j].chrom_idx;
                cur.first_snp = j;
            }
        }
    }
    return ranges;
}


// ===========================================================================
// Simplified consecutive state for permutation scans.
// We only track what is needed to decide run membership and update snp_freq.
// ===========================================================================

struct ConsecStatePerm {
    int32_t n_snp;
    int32_t n_opposite;
    int32_t n_missing;
    int32_t run_start_bp;
    int32_t run_end_bp;
    int32_t run_start_snp;  // local index (0-based within chromosome)
    int32_t run_end_snp;    // local index
    int32_t last_bp;
    uint8_t in_run;
};

static ConsecStatePerm make_consec_state_perm() {
    ConsecStatePerm s;
    memset(&s, 0, sizeof(s));
    return s;
}

// Emit a run: increment snp_freq_out for every SNP in [run_start_snp, run_end_snp].
static void emit_consec_perm(
    const ConsecStatePerm&  s,
    const ScanParams&       p,
    std::vector<int32_t>&   snp_freq_out)
{
    if (s.n_snp < p.min_snps) return;
    if ((s.run_end_bp - s.run_start_bp) < p.min_length_bp) return;
    for (int j = s.run_start_snp; j <= s.run_end_snp; ++j)
        snp_freq_out[j]++;
}

// Process one SNP (local index s within the chromosome) for one individual.
static void update_consec_perm(
    ConsecStatePerm&        st,
    int8_t                  geno,
    int32_t                 bp,
    int                     s,          // local SNP index
    const ScanParams&       p,
    std::vector<int32_t>&   snp_freq_out)
{
    // Gap check — no chromosome boundary here (single-chrom function)
    if (st.in_run && bp - st.last_bp >= p.max_gap) {
        emit_consec_perm(st, p, snp_freq_out);
        st.in_run = 0;
    }

    const int8_t opp = static_cast<int8_t>(1 - p.target);

    if (geno == p.target) {
        if (!st.in_run) {
            st.in_run        = 1;
            st.n_snp         = 0;
            st.n_opposite    = 0;
            st.n_missing     = 0;
            st.run_start_bp  = bp;
            st.run_start_snp = s;
            st.run_end_snp   = s;
        }
        st.n_snp++;
        st.run_end_bp  = bp;
        st.run_end_snp = s;

    } else if (geno == opp) {
        if (st.in_run) {
            if (st.n_opposite < p.max_opposite) {
                st.n_opposite++;
                st.n_snp++;
                st.run_end_bp  = bp;
                st.run_end_snp = s;
            } else {
                emit_consec_perm(st, p, snp_freq_out);
                st.in_run = 0;
            }
        }
    } else {
        // GENO_MISSING
        if (st.in_run) {
            if (st.n_missing < p.max_missing) {
                st.n_missing++;
                st.n_snp++;
                st.run_end_bp  = bp;
                st.run_end_snp = s;
            } else {
                emit_consec_perm(st, p, snp_freq_out);
                st.in_run = 0;
            }
        }
    }

    st.last_bp = bp;
}

// Scan one chromosome with the consecutive method and a permuted sample mapping.
// perm[i] = which BED-decoded slot individual i reads from.
// Fills snp_freq_out[0..n_c-1] (caller must zero-initialise before calling).
static void scan_chrom_consec_perm(
    const BedFile&          bed,
    const BimData&          bim,
    int                     first_snp,
    int                     n_c,
    int                     N,
    const std::vector<int>& perm,
    const ScanParams&       params,
    std::vector<int32_t>&   snp_freq_out)
{
    std::vector<ConsecStatePerm> states(N);
    for (int i = 0; i < N; ++i) states[i] = make_consec_state_perm();

    std::vector<int8_t> row_buf(N);

    for (int s = 0; s < n_c; ++s) {
        decode_snp_row(bed, first_snp + s, row_buf.data());
        const int32_t bp = bim.snps[first_snp + s].bp_pos;
        for (int i = 0; i < N; ++i)
            update_consec_perm(states[i], row_buf[perm[i]], bp, s, params, snp_freq_out);
    }

    // Flush open runs at chromosome end
    for (int i = 0; i < N; ++i)
        if (states[i].in_run)
            emit_consec_perm(states[i], params, snp_freq_out);
}


// ===========================================================================
// Simplified sliding-window state for permutation scans.
// ===========================================================================

// Decide whether SNP at local index s is in ROH and update snp_freq_out.
// Uses the same ring-buffer logic as decide_snp_sw in scan_roh.cpp but
// only increments snp_freq rather than emitting a full RohRecord.

struct SlidingRunPerm {
    int     run_start_snp;  // local index of first SNP in current run
    int32_t run_start_bp;
    int32_t run_end_bp;
    int     n_snp_run;
    uint8_t in_run;
};

static void decide_and_update_sw_perm(
    IndivStateSW&           st,
    SlidingRunPerm&         run,
    int                     s,          // local SNP index being decided
    int                     from_win,
    int                     to_win,
    const ScanParams&       p,
    std::vector<int32_t>&   snp_freq_out)
{
    int h = 0;
    for (int w = from_win; w <= to_win; ++w)
        h += st.win_pass_buf[w % (MAX_WINDOW * 2)];

    const int  n_cov  = to_win - from_win + 1;
    const bool in_roh = (static_cast<float>(h) / n_cov > p.threshold);

    if (in_roh) {
        if (!run.in_run) {
            run.in_run        = 1;
            run.run_start_snp = s;
            run.run_start_bp  = st.bp_buf[s % MAX_WINDOW];
            run.run_end_bp    = run.run_start_bp;
            run.n_snp_run     = 0;
        }
        run.run_end_bp = st.bp_buf[s % MAX_WINDOW];
        run.n_snp_run++;
        // tentatively mark; only finalised on run close
    } else {
        if (run.in_run) {
            // Close run: apply min_snps / min_length_bp and write snp_freq
            if (run.n_snp_run >= p.min_snps &&
                (run.run_end_bp - run.run_start_bp) >= p.min_length_bp) {
                for (int k = run.run_start_snp; k < run.run_start_snp + run.n_snp_run; ++k)
                    snp_freq_out[k]++;
            }
            run.in_run = 0;
        }
    }
}

static void scan_chrom_sliding_perm(
    const BedFile&          bed,
    const BimData&          bim,
    int                     first_snp,
    int                     n_c,
    int                     N,
    const std::vector<int>& perm,
    const ScanParams&       params,
    std::vector<int32_t>&   snp_freq_out)
{
    if (params.window_size > MAX_WINDOW)
        throw std::runtime_error("window_size exceeds MAX_WINDOW (256)");
    if (n_c == 0) return;

    const int    W   = params.window_size;
    const int8_t opp = static_cast<int8_t>(1 - params.target);

    std::vector<IndivStateSW>    states(N);
    std::vector<SlidingRunPerm>  runs(N);
    for (int i = 0; i < N; ++i) {
        memset(&states[i], 0, sizeof(IndivStateSW));
        memset(&runs[i],   0, sizeof(SlidingRunPerm));
    }

    std::vector<int8_t> row_buf(N);

    for (int s = 0; s < n_c; ++s) {
        decode_snp_row(bed, first_snp + s, row_buf.data());
        const int32_t bp = bim.snps[first_snp + s].bp_pos;

        for (int i = 0; i < N; ++i) {
            IndivStateSW&   st  = states[i];
            const int       pos = st.chrom_pos;
            const int8_t    g   = row_buf[perm[i]];

            st.geno_buf[pos % MAX_WINDOW] = g;
            st.bp_buf  [pos % MAX_WINDOW] = bp;

            if (pos >= W - 1) {
                const int w = pos - W + 1;

                // Evaluate window w = [w .. pos]
                bool gap_fail = false;
                for (int k = w; k < pos; ++k) {
                    if (st.bp_buf[(k + 1) % MAX_WINDOW] -
                        st.bp_buf[k % MAX_WINDOW] > params.max_gap) {
                        gap_fail = true;
                        break;
                    }
                }

                int8_t win_pass = 0;
                if (!gap_fail) {
                    int n_opp = 0, n_miss = 0;
                    for (int k = w; k <= pos; ++k) {
                        const int8_t gk = st.geno_buf[k % MAX_WINDOW];
                        if      (gk == GENO_MISSING) ++n_miss;
                        else if (gk == opp)          ++n_opp;
                    }
                    win_pass = (n_opp <= params.max_opposite &&
                                n_miss <= params.max_missing) ? 1 : 0;
                }
                st.win_pass_buf[w % (MAX_WINDOW * 2)] = win_pass;

                const int from_w = (w - W + 1 > 0) ? (w - W + 1) : 0;
                decide_and_update_sw_perm(st, runs[i], w, from_w, w,
                                          params, snp_freq_out);
            }
            st.chrom_pos++;
        }
    }

    // Flush trailing SNPs (last W-1 undecided)
    for (int i = 0; i < N; ++i) {
        IndivStateSW&  st  = states[i];
        SlidingRunPerm& run = runs[i];
        const int C = st.chrom_pos;
        if (C >= W) {
            const int last_win = C - W;
            for (int s = last_win + 1; s < C; ++s) {
                const int from_w = (s - W + 1 > 0) ? (s - W + 1) : 0;
                decide_and_update_sw_perm(st, run, s, from_w, last_win,
                                          params, snp_freq_out);
            }
        }
        // Close any still-open run
        if (run.in_run) {
            if (run.n_snp_run >= params.min_snps &&
                (run.run_end_bp - run.run_start_bp) >= params.min_length_bp) {
                for (int k = run.run_start_snp; k < run.run_start_snp + run.n_snp_run; ++k)
                    snp_freq_out[k]++;
            }
            run.in_run = 0;
        }
    }
}


// ===========================================================================
// Main entry point
// ===========================================================================

PermResult permutation_roh_islands(
    const BedFile&              bed,
    const BimData&              bim,
    const FamData&              fam,
    const std::vector<int32_t>& real_freq,
    const ScanParams&           params,
    int                         n_perm,
    double                      percentile,
    uint32_t                    seed)
{
    const int N = static_cast<int>(fam.samples.size());
    const int M = static_cast<int>(bim.snps.size());

    auto chrom_ranges = build_chrom_ranges(bim);
    const int n_chroms = static_cast<int>(chrom_ranges.size());

    PermResult result;
    result.chrom_indices.resize(n_chroms);
    result.thresholds.resize(n_chroms, 0.0);
    result.is_island.assign(M, false);

    if (N == 0 || M == 0 || n_perm <= 0) return result;

    // --- Generate all permutations upfront (MT19937, reproducible) ---
    std::mt19937 master_rng(seed != 0u ? seed : std::random_device{}());
    std::vector<std::vector<int>> perms(n_perm, std::vector<int>(N));
    {
        std::vector<int> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        for (int p = 0; p < n_perm; ++p) {
            perms[p] = idx;
            std::shuffle(perms[p].begin(), perms[p].end(), master_rng);
        }
    }

#ifdef _OPENMP
    const int actual = (params.n_threads > 0) ? params.n_threads : 1;
#else
    const int actual = 1;
    (void)actual;   // suppress unused-variable warning when OpenMP is absent
#endif

    // --- Process one chromosome at a time to bound peak memory ---
    for (int c = 0; c < n_chroms; ++c) {
        const ChromRange& cr = chrom_ranges[c];
        result.chrom_indices[c] = cr.chrom_idx;

        // pool[s * n_perm + p] = SNPROH count for local SNP s in permutation p
        // int16_t: SNPROH <= N samples < 32767 for any realistic dataset
        const size_t pool_size = static_cast<size_t>(cr.n_snps) * n_perm;
        std::vector<int16_t> pool(pool_size, 0);

        Rprintf("  Chromosome %d/%d (index %u, %d SNPs): running %d permutations ...",
                c + 1, n_chroms, static_cast<unsigned>(cr.chrom_idx),
                cr.n_snps, n_perm);
        R_FlushConsole();

#ifdef _OPENMP
        #pragma omp parallel for num_threads(actual) schedule(dynamic, 4)
#endif
        for (int p = 0; p < n_perm; ++p) {
            std::vector<int32_t> freq_local(cr.n_snps, 0);

            if (params.method == 0)
                scan_chrom_consec_perm(bed, bim, cr.first_snp, cr.n_snps, N,
                                       perms[p], params, freq_local);
            else
                scan_chrom_sliding_perm(bed, bim, cr.first_snp, cr.n_snps, N,
                                        perms[p], params, freq_local);

            for (int s = 0; s < cr.n_snps; ++s)
                pool[static_cast<size_t>(s) * n_perm + p] =
                    static_cast<int16_t>(freq_local[s]);
        }

        // Compute percentile of the full pool (sort in-place)
        std::sort(pool.begin(), pool.end());
        size_t idx_pct = static_cast<size_t>(percentile * static_cast<double>(pool_size));
        if (idx_pct >= pool_size) idx_pct = pool_size - 1u;
        result.thresholds[c] = static_cast<double>(pool[idx_pct]);

        // Flag real SNPs above threshold
        for (int s = 0; s < cr.n_snps; ++s) {
            int j = cr.first_snp + s;
            if (static_cast<double>(real_freq[j]) > result.thresholds[c])
                result.is_island[j] = true;
        }

        Rprintf(" done (threshold = %.1f)\n", result.thresholds[c]);
        R_FlushConsole();
    }

    return result;
}
