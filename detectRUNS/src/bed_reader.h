#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

#ifdef _WIN32
#  include <windows.h>
#endif


// ---------------------------------------------------------------------------
// Genotype values used throughout the scan engine
// ---------------------------------------------------------------------------

static const int8_t GENO_HOM     =  0;
static const int8_t GENO_HET     =  1;
static const int8_t GENO_MISSING = -1;


// ---------------------------------------------------------------------------
// BED layout flag
// ---------------------------------------------------------------------------

enum class BedLayout : uint8_t {
    SNP_MAJOR        = 0x01,
    INDIVIDUAL_MAJOR = 0x00
};


// ---------------------------------------------------------------------------
// SNP metadata (one entry per row in the BIM file)
// ---------------------------------------------------------------------------

struct SnpInfo {
    uint8_t     chrom_idx;   // mapped from chromosome name string
    int32_t     bp_pos;      // base-pair position
    std::string snp_name;    // SNP identifier (not used in hot loop)
};


// ---------------------------------------------------------------------------
// Sample metadata (one entry per row in the FAM file)
// ---------------------------------------------------------------------------

struct SampleInfo {
    std::string fid;
    std::string iid;
};


// ---------------------------------------------------------------------------
// Parsed BIM file
// ---------------------------------------------------------------------------

struct BimData {
    std::vector<SnpInfo>                        snps;
    std::unordered_map<std::string, uint8_t>    chrom_map;    // name → index
    std::vector<std::string>                    chrom_names;  // index → name
};


// ---------------------------------------------------------------------------
// Parsed FAM file
// ---------------------------------------------------------------------------

struct FamData {
    std::vector<SampleInfo> samples;
};


// ---------------------------------------------------------------------------
// Open BED file (mmap-backed)
// ---------------------------------------------------------------------------

struct BedFile {
    const uint8_t* data      = nullptr;
    size_t         file_size = 0;
    BedLayout      layout    = BedLayout::SNP_MAJOR;
    int            n_samples = 0;
    int            n_snps    = 0;
    int            row_bytes = 0;   // bytes per SNP row   (SNP-major)
    int            col_bytes = 0;   // bytes per indiv row (individual-major)

#ifdef _WIN32
    HANDLE win_file = INVALID_HANDLE_VALUE;
    HANDLE win_map  = nullptr;
#else
    int fd = -1;
#endif
};


// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

BimData parse_bim(const std::string& bim_path);
FamData parse_fam(const std::string& fam_path);

BedFile open_bed(const std::string& bed_path, int n_samples, int n_snps);
void    close_bed(BedFile& bed);

// Decode one full SNP row (SNP-major) into out[0..n_samples-1]
void decode_snp_row(const BedFile& bed, int snp_idx, int8_t* out);

// Decode one full individual row (individual-major) into out[0..n_snps-1]
void decode_individual_row(const BedFile& bed, int ind_idx, int8_t* out);
