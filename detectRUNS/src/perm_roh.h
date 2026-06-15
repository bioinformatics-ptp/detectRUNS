#pragma once

#include "scan_roh.h"
#include <vector>
#include <cstdint>


// ---------------------------------------------------------------------------
// ChromRange: contiguous block of SNPs belonging to one chromosome.
// Indices are 0-based positions into BimData::snps.
// ---------------------------------------------------------------------------

struct ChromRange {
    uint8_t chrom_idx;
    int     first_snp;   // first SNP index in this chromosome (inclusive)
    int     n_snps;      // number of SNPs on this chromosome
};


// ---------------------------------------------------------------------------
// PermResult: output of permutation_roh_islands()
//
// chrom_indices[k] / thresholds[k] are parallel vectors: the internal
// uint8_t chromosome index and its derived SNPROH threshold.
// is_island[j]  (size == n_snps) is true if real snp_freq[j] > threshold
// for the chromosome j belongs to.
// ---------------------------------------------------------------------------

struct PermResult {
    std::vector<uint8_t> chrom_indices;
    std::vector<double>  thresholds;
    std::vector<bool>    is_island;
};


// ---------------------------------------------------------------------------
// Main permutation entry point.
//
// For each chromosome independently:
//   1. Generate n_perm random permutations of sample indices.
//   2. Run the ROH scan (same method / params as the original scan) on each
//      permuted dataset, recording per-SNP SNPROH counts.
//   3. Pool all SNP x perm SNPROH values for that chromosome.
//   4. Take the 'percentile'-th quantile as the chromosome-specific threshold.
//
// Parallelism (OpenMP, when compiled with -fopenmp): permutations are
// distributed across threads; each thread works on its own private buffers.
//
// Parameters
//   bed         open BedFile (mmap)
//   bim         parsed BIM metadata
//   fam         parsed FAM metadata
//   real_freq   observed SNPROH counts from the original scan (size n_snps)
//   params      identical ScanParams used for the original scan
//   n_perm      number of permutations
//   percentile  quantile for threshold (e.g. 0.99 for the 99th percentile)
//   seed        MT19937 seed; 0 = draw from std::random_device
// ---------------------------------------------------------------------------

PermResult permutation_roh_islands(
    const BedFile&              bed,
    const BimData&              bim,
    const FamData&              fam,
    const std::vector<int32_t>& real_freq,
    const ScanParams&           params,
    int                         n_perm,
    double                      percentile,
    uint32_t                    seed
);
