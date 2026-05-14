#pragma once

#include "scan_roh.h"
#include <string>
#include <vector>


// ---------------------------------------------------------------------------
// Data returned by read_roh_binary — fully self-contained
// ---------------------------------------------------------------------------

struct RohBinaryData {
    std::vector<RohRecord>   records;
    FamData                  fam;
    std::vector<std::string> chrom_names;  // sparse; indexed by chrom_idx (0..255)
};


// ---------------------------------------------------------------------------
// Binary file I/O
//
// write_roh_binary: serialises records + sample names + chromosome names.
//   File magic "ROHB", version 1, little-endian integers, packed RohRecord.
//
// read_roh_binary: reads a file written by write_roh_binary and returns
//   a RohBinaryData with records, fam, and chrom_names fully reconstructed.
//
// Both functions throw std::runtime_error on I/O or format errors.
// ---------------------------------------------------------------------------

void write_roh_binary(
    const std::string&            path,
    const std::vector<RohRecord>& records,
    const FamData&                fam,
    const BimData&                bim
);

RohBinaryData read_roh_binary(const std::string& path);
