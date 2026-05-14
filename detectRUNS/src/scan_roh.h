#pragma once

#include "bed_reader.h"
#include <vector>
#include <cstdint>


// ---------------------------------------------------------------------------
// ROH record: one detected run, per individual (21 bytes, packed)
// ---------------------------------------------------------------------------

struct RohRecord {
    int32_t  sample_idx;   // index into FamData::samples
    uint8_t  chrom_idx;    // chromosome index (from BimData::chrom_map)
    int32_t  start_bp;     // start position in base pairs
    int32_t  end_bp;       // end position in base pairs
    int32_t  n_snps;       // SNPs inside run (target + accepted opp + accepted miss)
    int16_t  n_het;        // opposite-type genotypes inside run
    int16_t  n_missing;    // missing genotypes inside run
} __attribute__((packed));

static_assert(sizeof(RohRecord) == 21, "RohRecord must be exactly 21 bytes");


// ---------------------------------------------------------------------------
// Per-individual summary accumulated during the scan
// ---------------------------------------------------------------------------

struct IndivSummary {
    int32_t  n_roh;
    int64_t  total_length_bp;
    int32_t  n_snps_in_roh;
};


// ---------------------------------------------------------------------------
// All scan parameters in one struct
// ---------------------------------------------------------------------------

struct ScanParams {
    int8_t   target;        // 0 = ROHom (GENO_HOM), 1 = ROHet (GENO_HET)
    int      method;        // 0 = consecutive, 1 = sliding window
    int      min_snps;      // minimum SNPs in a qualifying run
    int      max_opposite;  // consecutive: max opp SNPs in run; sliding: max opp in window
    int      max_missing;   // consecutive: max missing in run;  sliding: max missing in window
    int      min_length_bp; // minimum run length in base pairs
    int      max_gap;       // max gap between consecutive SNPs (bp); >= breaks a run
    int      window_size;   // sliding window only: window width in SNPs
    double   threshold;     // sliding window only: Bjelland coverage ratio threshold (>)
    int      n_threads;     // number of OpenMP threads (1 = single-threaded)
    bool     verbose;       // print progress bar and per-chromosome updates
};


// ---------------------------------------------------------------------------
// Streaming sliding-window state (Phase 5) — one per individual
//
// Ring buffers keep only the last MAX_WINDOW SNPs and 2×MAX_WINDOW window
// results; O(N × W) memory rather than O(N × max_chrom_snps).
// ---------------------------------------------------------------------------

#define MAX_WINDOW 256

struct IndivStateSW {
    int8_t   geno_buf[MAX_WINDOW];        // ring: genotypes (last W SNPs, within-chrom)
    int32_t  bp_buf[MAX_WINDOW];          // ring: bp positions (last W SNPs)
    int8_t   win_pass_buf[MAX_WINDOW * 2];// ring: window pass/fail (last 2W windows)
    int32_t  chrom_pos;                   // SNP index within current chromosome (0-based)
    int32_t  n_opp_run;                   // opposite genotypes accumulated in current run
    int32_t  n_miss_run;                  // missing genotypes accumulated in current run
    int32_t  n_snp_run;                   // total SNPs in current run
    int32_t  run_start_bp;
    int32_t  run_end_bp;
    int32_t  run_start_snp;               // global BIM index of first SNP in run
    uint8_t  in_run;
};


// ---------------------------------------------------------------------------
// Main scan entry point
//
// Reads the BED file (SNP-major layout) and detects ROH for all individuals.
// Populates:
//   records   - one RohRecord per detected run
//   summaries - one IndivSummary per individual (indexed by FAM order)
//   snp_freq  - one counter per SNP: times the SNP fell inside any ROH
//
// Throws std::runtime_error on unsupported BED layout or invalid parameters.
// ---------------------------------------------------------------------------

void scan_roh_bed(
    const BedFile&             bed,
    const BimData&             bim,
    const FamData&             fam,
    const ScanParams&          params,
    std::vector<RohRecord>&    records,
    std::vector<IndivSummary>& summaries,
    std::vector<int32_t>&      snp_freq
);
