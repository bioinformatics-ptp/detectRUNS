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
    const int8_t*           chrom_genos,   // [n_c × N] row-major, pre-decoded
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

    for (int s = 0; s < n_c; ++s) {
        const int8_t* row = chrom_genos + static_cast<size_t>(s) * N;
        const int32_t bp  = bim.snps[first_snp + s].bp_pos;
        for (int i = 0; i < N; ++i)
            update_consec_perm(states[i], row[perm[i]], bp, s, params, snp_freq_out);
    }

    // Flush open runs at chromosome end
    for (int i = 0; i < N; ++i)
        if (states[i].in_run)
            emit_consec_perm(states[i], params, snp_freq_out);
}


// ===========================================================================
// Optimised sliding-window permutation scan.
//
// The original implementation re-evaluated each window from scratch with
// three O(windowSize) inner loops:
//   1. gap check             — O(W) scan
//   2. n_opp / n_miss count  — O(W) scan
//   3. win_pass_buf sum      — O(W) scan
//
// This version replaces all three with O(1) sliding-sum updates:
//   • n_opp_win / n_miss_win: add incoming SNP, drop departing SNP
//   • n_bad_gaps:             add new gap, drop departing gap
//   • win_pass_sum:           add new win_pass, subtract the one that
//                             left the per-SNP covering window set
//
// Expected wall-clock speedup ≈ window_size (≈ 20× for default W=20).
// The trailing W-1 SNPs still use an O(W^2) loop, but that is at most
// W-1 iterations total per individual so its cost is negligible.
// ===========================================================================

static void scan_chrom_sliding_perm(
    const int8_t*           chrom_genos,   // [n_c × N] row-major, pre-decoded
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

    // Per-individual state — all sliding sums live here.
    struct IndivFast {
        int8_t  geno_buf    [MAX_WINDOW * 2];   // ring: genotypes
        int32_t bp_buf      [MAX_WINDOW * 2];   // ring: base-pair positions
        int8_t  win_pass_buf[MAX_WINDOW * 2];   // ring: per-window pass flag
        int     n_opp_win;      // opposite genotypes inside current window
        int     n_miss_win;     // missing   genotypes inside current window
        int     n_bad_gaps;     // gaps > max_gap    inside current window
        int     win_pass_sum;   // sliding sum of win_pass over the W windows
                                // that cover the SNP currently being decided
        int     chrom_pos;
        // run tracking
        int     run_start_snp;
        int32_t run_start_bp;
        int32_t run_end_bp;
        int     n_snp_run;
        bool    in_run;
    };

    std::vector<IndivFast> states(N);
    for (auto& st : states) memset(&st, 0, sizeof(IndivFast));

    for (int s = 0; s < n_c; ++s) {
        const int8_t* row = chrom_genos + static_cast<size_t>(s) * N;
        const int32_t bp  = bim.snps[first_snp + s].bp_pos;

        for (int i = 0; i < N; ++i) {
            IndivFast&   st  = states[i];
            const int    pos = st.chrom_pos;
            const int8_t g   = row[perm[i]];
            const int    ri  = pos % (MAX_WINDOW * 2);   // ring index

            st.geno_buf[ri] = g;
            st.bp_buf  [ri] = bp;

            // --- O(1): add incoming SNP to window counts ---
            if      (g == opp)          ++st.n_opp_win;
            else if (g == GENO_MISSING) ++st.n_miss_win;

            // new gap: (pos-1) → pos
            if (pos > 0) {
                const int32_t prev_bp = st.bp_buf[(pos - 1) % (MAX_WINDOW * 2)];
                if (bp - prev_bp > params.max_gap) ++st.n_bad_gaps;
            }

            // --- O(1): drop departing SNP once the window is full ---
            if (pos >= W) {
                const int dep   = pos - W;
                const int8_t dg = st.geno_buf[dep % (MAX_WINDOW * 2)];
                if      (dg == opp)          --st.n_opp_win;
                else if (dg == GENO_MISSING) --st.n_miss_win;

                // departing gap: dep → dep+1
                const int32_t dep_bp = st.bp_buf[ dep      % (MAX_WINDOW * 2)];
                const int32_t nxt_bp = st.bp_buf[(dep + 1) % (MAX_WINDOW * 2)];
                if (nxt_bp - dep_bp > params.max_gap) --st.n_bad_gaps;
            }

            // --- Evaluate window once it is complete ---
            if (pos >= W - 1) {
                const int w = pos - W + 1;   // window-start SNP index

                // O(1) window pass flag
                const int8_t win_pass = (st.n_bad_gaps   == 0 &&
                                         st.n_opp_win  <= params.max_opposite &&
                                         st.n_miss_win <= params.max_missing) ? 1 : 0;
                st.win_pass_buf[w % (MAX_WINDOW * 2)] = win_pass;

                // O(1) sliding sum of covering windows for SNP at position w
                st.win_pass_sum += win_pass;
                if (w >= W)
                    st.win_pass_sum -= st.win_pass_buf[(w - W) % (MAX_WINDOW * 2)];

                // Decide whether SNP w is in ROH
                const int  n_cov  = (w < W) ? (w + 1) : W;
                const bool in_roh = (static_cast<float>(st.win_pass_sum) / n_cov
                                     > params.threshold);

                if (in_roh) {
                    if (!st.in_run) {
                        st.in_run        = true;
                        st.run_start_snp = w;
                        st.run_start_bp  = st.bp_buf[w % (MAX_WINDOW * 2)];
                        st.run_end_bp    = st.run_start_bp;
                        st.n_snp_run     = 0;
                    }
                    st.run_end_bp = st.bp_buf[w % (MAX_WINDOW * 2)];
                    st.n_snp_run++;
                } else {
                    if (st.in_run) {
                        if (st.n_snp_run >= params.min_snps &&
                            (st.run_end_bp - st.run_start_bp) >= params.min_length_bp)
                            for (int k = st.run_start_snp;
                                     k < st.run_start_snp + st.n_snp_run; ++k)
                                snp_freq_out[k]++;
                        st.in_run = false;
                    }
                }
            }

            st.chrom_pos++;
        }
    }

    // Flush trailing W-1 undecided SNPs (at most W-1 per individual — negligible cost)
    for (int i = 0; i < N; ++i) {
        IndivFast& st = states[i];
        const int  C  = st.chrom_pos;

        if (C < W) {
            // Chromosome shorter than one window — nothing was ever decided.
            if (st.in_run) {
                if (st.n_snp_run >= params.min_snps &&
                    (st.run_end_bp - st.run_start_bp) >= params.min_length_bp)
                    for (int k = st.run_start_snp;
                             k < st.run_start_snp + st.n_snp_run; ++k)
                        snp_freq_out[k]++;
            }
            continue;
        }

        const int last_win = C - W;   // last window evaluated in the main loop

        // SNPs last_win+1 … C-1 (at most W-1 of them): use simple O(W) loop —
        // total iterations <= (W-1)^2 / 2, completely dominated by the main loop.
        for (int q = last_win + 1; q < C; ++q) {
            const int from_w = (q - W + 1 > 0) ? (q - W + 1) : 0;
            int h = 0;
            for (int ww = from_w; ww <= last_win; ++ww)
                h += st.win_pass_buf[ww % (MAX_WINDOW * 2)];
            const int n_cov = last_win - from_w + 1;
            if (n_cov <= 0) continue;
            const bool in_roh = (static_cast<float>(h) / n_cov > params.threshold);
            if (in_roh) {
                if (!st.in_run) {
                    st.in_run        = true;
                    st.run_start_snp = q;
                    st.run_start_bp  = st.bp_buf[q % (MAX_WINDOW * 2)];
                    st.run_end_bp    = st.run_start_bp;
                    st.n_snp_run     = 0;
                }
                st.run_end_bp = st.bp_buf[q % (MAX_WINDOW * 2)];
                st.n_snp_run++;
            } else {
                if (st.in_run) {
                    if (st.n_snp_run >= params.min_snps &&
                        (st.run_end_bp - st.run_start_bp) >= params.min_length_bp)
                        for (int k = st.run_start_snp;
                                 k < st.run_start_snp + st.n_snp_run; ++k)
                            snp_freq_out[k]++;
                    st.in_run = false;
                }
            }
        }

        if (st.in_run) {
            if (st.n_snp_run >= params.min_snps &&
                (st.run_end_bp - st.run_start_bp) >= params.min_length_bp)
                for (int k = st.run_start_snp;
                         k < st.run_start_snp + st.n_snp_run; ++k)
                    snp_freq_out[k]++;
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

        // Decode every SNP row for this chromosome once; permutations reuse the cache.
        std::vector<int8_t> chrom_genos(static_cast<size_t>(cr.n_snps) * N);
        for (int s = 0; s < cr.n_snps; ++s)
            decode_snp_row(bed, cr.first_snp + s,
                           chrom_genos.data() + static_cast<size_t>(s) * N);

#ifdef _OPENMP
        #pragma omp parallel for num_threads(actual) schedule(dynamic, 4)
#endif
        for (int p = 0; p < n_perm; ++p) {
            std::vector<int32_t> freq_local(cr.n_snps, 0);

            if (params.method == 0)
                scan_chrom_consec_perm(chrom_genos.data(), bim, cr.first_snp, cr.n_snps, N,
                                       perms[p], params, freq_local);
            else
                scan_chrom_sliding_perm(chrom_genos.data(), bim, cr.first_snp, cr.n_snps, N,
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
