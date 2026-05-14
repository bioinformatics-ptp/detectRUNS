#include "bed_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#endif


// ---------------------------------------------------------------------------
// 256-entry LUT: maps one BED byte to four int8_t genotype values.
//
// BED 2-bit encoding (two bits per sample, LSB first within byte):
//   00 = homozygous first allele  → 0
//   01 = missing                  → GENO_MISSING (-1)
//   10 = heterozygous             → 1
//   11 = homozygous second allele → 0
//
// The LUT is indexed by the raw byte value (0-255). Each entry is an array
// of four genotype values for the four samples packed in that byte.
// ---------------------------------------------------------------------------

static bool lut_ready = false;
static int8_t BED_LUT[256][4];

static void build_lut() {
    if (lut_ready) return;

    static const int8_t decode[4] = {
        GENO_HOM,      // 00 → homozygous
        GENO_MISSING,  // 01 → missing
        GENO_HET,      // 10 → heterozygous
        GENO_HOM       // 11 → homozygous
    };

    for (int byte_val = 0; byte_val < 256; byte_val++) {
        for (int sample_in_byte = 0; sample_in_byte < 4; sample_in_byte++) {
            int two_bits = (byte_val >> (sample_in_byte * 2)) & 0x03;
            BED_LUT[byte_val][sample_in_byte] = decode[two_bits];
        }
    }
    lut_ready = true;
}


// ---------------------------------------------------------------------------
// Chromosome name → index mapping
// Numeric chromosomes map to their integer value (1-based).
// Named chromosomes get indices starting at CHROM_NAMED_START.
// ---------------------------------------------------------------------------

static const uint8_t CHROM_NAMED_START = 100;

static const struct { const char* name; uint8_t idx; } NAMED_CHROMS[] = {
    {"X",  100}, {"Y",  101}, {"XY", 102},
    {"MT", 103}, {"M",  103},
    {"Z",  104}, {"W",  105},
    {nullptr, 0}
};

static uint8_t chrom_name_to_idx(const std::string& name,
                                  std::unordered_map<std::string, uint8_t>& map,
                                  std::vector<std::string>& rev) {
    auto it = map.find(name);
    if (it != map.end()) return it->second;

    // try numeric first
    bool is_numeric = !name.empty() &&
        std::all_of(name.begin(), name.end(), ::isdigit);
    if (is_numeric) {
        int n = std::stoi(name);
        uint8_t idx = static_cast<uint8_t>(n);
        map[name] = idx;
        if (idx >= rev.size()) rev.resize(idx + 1);
        rev[idx] = name;
        return idx;
    }

    // try named table
    for (int i = 0; NAMED_CHROMS[i].name != nullptr; i++) {
        if (name == NAMED_CHROMS[i].name) {
            uint8_t idx = NAMED_CHROMS[i].idx;
            map[name] = idx;
            if (idx >= rev.size()) rev.resize(idx + 1);
            rev[idx] = name;
            return idx;
        }
    }

    // unknown — assign next available index above CHROM_NAMED_START
    uint8_t idx = static_cast<uint8_t>(CHROM_NAMED_START + rev.size());
    map[name] = idx;
    rev.push_back(name);
    return idx;
}


// ---------------------------------------------------------------------------
// BIM parser
// Format: CHR  SNP_ID  cM  BP  A1  A2  (tab or space separated)
// ---------------------------------------------------------------------------

BimData parse_bim(const std::string& bim_path) {
    BimData bim;
    std::ifstream f(bim_path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open BIM file: " + bim_path);

    std::string line, chrom_str, snp_id, cm_str, bp_str, a1, a2;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        if (!(ss >> chrom_str >> snp_id >> cm_str >> bp_str)) continue;
        ss >> a1 >> a2;  // optional — not used in scan

        SnpInfo si;
        si.chrom_idx = chrom_name_to_idx(chrom_str, bim.chrom_map, bim.chrom_names);
        si.bp_pos    = std::stoi(bp_str);
        si.snp_name  = snp_id;
        bim.snps.push_back(si);
    }
    return bim;
}


// ---------------------------------------------------------------------------
// FAM parser
// Format: FID  IID  PAT  MAT  SEX  PHENO  (space separated)
// ---------------------------------------------------------------------------

FamData parse_fam(const std::string& fam_path) {
    FamData fam;
    std::ifstream f(fam_path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open FAM file: " + fam_path);

    std::string line, fid, iid, rest;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        if (!(ss >> fid >> iid)) continue;
        SampleInfo si;
        si.fid = fid;
        si.iid = iid;
        fam.samples.push_back(si);
    }
    return fam;
}


// ---------------------------------------------------------------------------
// BED file open: validate magic bytes, detect layout
// ---------------------------------------------------------------------------

BedFile open_bed(const std::string& bed_path,
                 int n_samples, int n_snps) {
    BedFile bed;
    bed.n_samples  = n_samples;
    bed.n_snps     = n_snps;
    bed.row_bytes  = (n_samples + 3) / 4;   // bytes per SNP row (SNP-major)
    bed.col_bytes  = (n_snps   + 3) / 4;    // bytes per individual row (individual-major)
    bed.file_size  = 3 + static_cast<size_t>(n_snps) * bed.row_bytes;

    std::ifstream f(bed_path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("Cannot open BED file: " + bed_path);

    uint8_t magic[3];
    f.read(reinterpret_cast<char*>(magic), 3);
    if (magic[0] != 0x6C || magic[1] != 0x1B)
        throw std::runtime_error("Not a valid PLINK BED file: " + bed_path);

    if (magic[2] == 0x01) {
        bed.layout = BedLayout::SNP_MAJOR;
    } else if (magic[2] == 0x00) {
        bed.layout = BedLayout::INDIVIDUAL_MAJOR;
    } else {
        throw std::runtime_error("Unknown BED layout byte in: " + bed_path);
    }
    f.close();

    // memory-map the file
#ifdef _WIN32
    bed.win_file = CreateFileA(bed_path.c_str(), GENERIC_READ, FILE_SHARE_READ,
                               nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (bed.win_file == INVALID_HANDLE_VALUE)
        throw std::runtime_error("CreateFile failed: " + bed_path);

    bed.win_map = CreateFileMappingA(bed.win_file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (!bed.win_map)
        throw std::runtime_error("CreateFileMapping failed: " + bed_path);

    bed.data = static_cast<const uint8_t*>(
        MapViewOfFile(bed.win_map, FILE_MAP_READ, 0, 0, 0));
    if (!bed.data)
        throw std::runtime_error("MapViewOfFile failed: " + bed_path);
#else
    bed.fd = open(bed_path.c_str(), O_RDONLY);
    if (bed.fd < 0)
        throw std::runtime_error("open() failed: " + bed_path);

    bed.data = static_cast<const uint8_t*>(
        mmap(nullptr, bed.file_size, PROT_READ, MAP_SHARED, bed.fd, 0));
    if (bed.data == MAP_FAILED)
        throw std::runtime_error("mmap() failed: " + bed_path);
#endif

    build_lut();
    return bed;
}


// ---------------------------------------------------------------------------
// BED file close: release mmap and handles
// ---------------------------------------------------------------------------

void close_bed(BedFile& bed) {
#ifdef _WIN32
    if (bed.data)     UnmapViewOfFile(bed.data);
    if (bed.win_map)  CloseHandle(bed.win_map);
    if (bed.win_file != INVALID_HANDLE_VALUE) CloseHandle(bed.win_file);
    bed.data     = nullptr;
    bed.win_map  = nullptr;
    bed.win_file = INVALID_HANDLE_VALUE;
#else
    if (bed.data && bed.data != MAP_FAILED)
        munmap(const_cast<uint8_t*>(bed.data), bed.file_size);
    if (bed.fd >= 0) close(bed.fd);
    bed.data = nullptr;
    bed.fd   = -1;
#endif
}


// ---------------------------------------------------------------------------
// Decode one SNP row (SNP-major layout) into caller-supplied buffer.
// snp_idx:    0-based SNP index
// out[]:      caller-allocated buffer of size n_samples
// ---------------------------------------------------------------------------

void decode_snp_row(const BedFile& bed, int snp_idx, int8_t* out) {
    const uint8_t* row = bed.data + 3 + static_cast<size_t>(snp_idx) * bed.row_bytes;

    int n = bed.n_samples;
    int full_bytes = n / 4;
    int remainder  = n % 4;

    for (int b = 0; b < full_bytes; b++) {
        const int8_t* decoded = BED_LUT[row[b]];
        out[b*4 + 0] = decoded[0];
        out[b*4 + 1] = decoded[1];
        out[b*4 + 2] = decoded[2];
        out[b*4 + 3] = decoded[3];
    }
    if (remainder) {
        const int8_t* decoded = BED_LUT[row[full_bytes]];
        for (int r = 0; r < remainder; r++)
            out[full_bytes*4 + r] = decoded[r];
    }
}


// ---------------------------------------------------------------------------
// Decode one individual row (individual-major layout) into caller buffer.
// ind_idx:    0-based individual index
// out[]:      caller-allocated buffer of size n_snps
// ---------------------------------------------------------------------------

void decode_individual_row(const BedFile& bed, int ind_idx, int8_t* out) {
    const uint8_t* row = bed.data + 3 + static_cast<size_t>(ind_idx) * bed.col_bytes;

    int n = bed.n_snps;
    int full_bytes = n / 4;
    int remainder  = n % 4;

    for (int b = 0; b < full_bytes; b++) {
        const int8_t* decoded = BED_LUT[row[b]];
        out[b*4 + 0] = decoded[0];
        out[b*4 + 1] = decoded[1];
        out[b*4 + 2] = decoded[2];
        out[b*4 + 3] = decoded[3];
    }
    if (remainder) {
        const int8_t* decoded = BED_LUT[row[full_bytes]];
        for (int r = 0; r < remainder; r++)
            out[full_bytes*4 + r] = decoded[r];
    }
}
