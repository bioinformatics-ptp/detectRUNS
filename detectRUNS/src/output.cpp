#include "output.h"
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <set>


// ===========================================================================
// Binary file format  (version 1)
//
//  HEADER      64 bytes, packed (see RohBinHeader below)
//  CHROM TABLE n_chroms entries:  uint8 idx | uint8 name_len | char name[]
//  SAMPLE TABLE n_samples entries: uint8 fid_len | char fid[] |
//                                  uint8 iid_len | char iid[]
//  ROH RECORDS n_records × 21 bytes  (packed RohRecord, see scan_roh.h)
//
//  All multi-byte integers are native (little-endian on x86 / arm64).
// ===========================================================================

#pragma pack(push, 1)
struct RohBinHeader {
    char     magic[4];      // "ROHB"
    uint32_t version;       // 1
    uint64_t n_records;
    uint32_t n_samples;
    uint32_t n_snps;        // BIM row count (for cross-validation)
    uint8_t  n_chroms;      // unique chromosomes in this file
    uint8_t  reserved[39];  // zeros
};
#pragma pack(pop)

static_assert(sizeof(RohBinHeader) == 64, "RohBinHeader must be 64 bytes");


// ---------------------------------------------------------------------------
// Low-level I/O helpers (throw on any failure)
// ---------------------------------------------------------------------------

static void w_raw(FILE* f, const void* buf, size_t n) {
    if (fwrite(buf, 1, n, f) != n)
        throw std::runtime_error("binary write error");
}
static void r_raw(FILE* f, void* buf, size_t n) {
    if (fread(buf, 1, n, f) != n)
        throw std::runtime_error("binary read error (unexpected end of file)");
}
static void w_u8(FILE* f, uint8_t v) { w_raw(f, &v, 1); }
static uint8_t r_u8(FILE* f) { uint8_t v; r_raw(f, &v, 1); return v; }

static void w_str(FILE* f, const std::string& s) {
    if (s.size() > 255)
        throw std::runtime_error("string too long for binary format (max 255 chars)");
    w_u8(f, static_cast<uint8_t>(s.size()));
    if (!s.empty()) w_raw(f, s.data(), s.size());
}
static std::string r_str(FILE* f) {
    uint8_t len = r_u8(f);
    if (len == 0) return {};
    std::string s(len, '\0');
    r_raw(f, &s[0], len);
    return s;
}


// ===========================================================================
// write_roh_binary
// ===========================================================================

void write_roh_binary(
    const std::string&            path,
    const std::vector<RohRecord>& records,
    const FamData&                fam,
    const BimData&                bim)
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot open for writing: " + path);

    try {
        // Collect unique chrom indices from records (sorted)
        std::set<uint8_t> seen;
        for (const auto& r : records) seen.insert(r.chrom_idx);
        std::vector<uint8_t> chrom_idxs(seen.begin(), seen.end());

        if (chrom_idxs.size() > 255)
            throw std::runtime_error("too many unique chromosomes for binary format");

        // Write header
        RohBinHeader hdr;
        memset(&hdr, 0, sizeof(hdr));
        hdr.magic[0] = 'R'; hdr.magic[1] = 'O';
        hdr.magic[2] = 'H'; hdr.magic[3] = 'B';
        hdr.version   = 1;
        hdr.n_records = static_cast<uint64_t>(records.size());
        hdr.n_samples = static_cast<uint32_t>(fam.samples.size());
        hdr.n_snps    = static_cast<uint32_t>(bim.snps.size());
        hdr.n_chroms  = static_cast<uint8_t>(chrom_idxs.size());
        w_raw(f, &hdr, sizeof(hdr));

        // Write chrom table
        for (uint8_t cidx : chrom_idxs) {
            const std::string& name =
                (static_cast<size_t>(cidx) < bim.chrom_names.size() &&
                 !bim.chrom_names[cidx].empty())
                ? bim.chrom_names[cidx]
                : std::to_string(cidx);
            w_u8(f, cidx);
            w_str(f, name);
        }

        // Write sample table
        for (const auto& s : fam.samples) {
            w_str(f, s.fid);
            w_str(f, s.iid);
        }

        // Write ROH records (bulk)
        if (!records.empty())
            w_raw(f, records.data(), sizeof(RohRecord) * records.size());

    } catch (...) {
        fclose(f);
        throw;
    }
    fclose(f);
}


// ===========================================================================
// read_roh_binary
// ===========================================================================

RohBinaryData read_roh_binary(const std::string& path)
{
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open for reading: " + path);

    RohBinaryData out;
    out.chrom_names.resize(256);  // indexed by uint8_t chrom_idx

    try {
        RohBinHeader hdr;
        r_raw(f, &hdr, sizeof(hdr));

        if (hdr.magic[0] != 'R' || hdr.magic[1] != 'O' ||
            hdr.magic[2] != 'H' || hdr.magic[3] != 'B')
            throw std::runtime_error("not a valid ROHB file: " + path);
        if (hdr.version != 1)
            throw std::runtime_error("unsupported ROHB version: " + path);

        // Read chrom table
        for (uint8_t i = 0; i < hdr.n_chroms; ++i) {
            uint8_t    cidx = r_u8(f);
            std::string name = r_str(f);
            out.chrom_names[cidx] = name;
        }

        // Read sample table
        out.fam.samples.resize(hdr.n_samples);
        for (uint32_t i = 0; i < hdr.n_samples; ++i) {
            out.fam.samples[i].fid = r_str(f);
            out.fam.samples[i].iid = r_str(f);
        }

        // Read ROH records (bulk)
        out.records.resize(static_cast<size_t>(hdr.n_records));
        if (hdr.n_records > 0)
            r_raw(f, out.records.data(),
                  sizeof(RohRecord) * static_cast<size_t>(hdr.n_records));

    } catch (...) {
        fclose(f);
        throw;
    }
    fclose(f);
    return out;
}
