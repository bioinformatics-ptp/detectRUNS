#include <Rcpp.h>
#include <unordered_map>
#include "scan_roh.h"
#include "output.h"
#include "perm_roh.h"


// ===========================================================================
// Output assembly helpers (Rcpp — R/C++ boundary only)
// ===========================================================================

// Build 'runs' data.frame: group | id | chrom | nSNP | from | to | lengthBps
// Column format matches the output of scanRUNS().
static Rcpp::DataFrame build_runs_df(
    const std::vector<RohRecord>& records,
    const FamData& fam,
    const BimData& bim)
{
    const int n = static_cast<int>(records.size());

    Rcpp::CharacterVector group(n), id(n), chrom(n);
    Rcpp::IntegerVector   nSNP(n), from(n), to(n), lengthBps(n), nHet(n), nMissing(n);

    for (int k = 0; k < n; ++k) {
        const RohRecord&  r = records[k];
        const SampleInfo& s = fam.samples[r.sample_idx];

        group[k] = s.fid;
        id[k]    = s.iid;

        chrom[k] = (static_cast<size_t>(r.chrom_idx) < bim.chrom_names.size() &&
                    !bim.chrom_names[static_cast<size_t>(r.chrom_idx)].empty())
                   ? bim.chrom_names[static_cast<size_t>(r.chrom_idx)]
                   : std::to_string(static_cast<int>(r.chrom_idx));

        nSNP[k]      = r.n_snps;
        from[k]      = r.start_bp;
        to[k]        = r.end_bp;
        lengthBps[k] = r.end_bp - r.start_bp;
        nHet[k]      = r.n_het;
        nMissing[k]  = r.n_missing;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("group")     = group,
        Rcpp::Named("id")        = id,
        Rcpp::Named("chrom")     = chrom,
        Rcpp::Named("nSNP")      = nSNP,
        Rcpp::Named("from")      = from,
        Rcpp::Named("to")        = to,
        Rcpp::Named("lengthBps") = lengthBps,
        Rcpp::Named("nHet")      = nHet,
        Rcpp::Named("nMissing")  = nMissing,
        Rcpp::_["stringsAsFactors"] = false
    );
}

// Build per-individual 'summary' data.frame.
// n_ROH, total_length_bp, mean_length, n_snps_in_roh are raw values;
// R wrapper adds F_ROH using the total autosomal genome length.
static Rcpp::DataFrame build_summary_df(
    const std::vector<IndivSummary>& summaries,
    const FamData& fam)
{
    const int N = static_cast<int>(summaries.size());

    Rcpp::CharacterVector group(N), id(N);
    Rcpp::IntegerVector   n_roh(N), n_snps_in_roh(N);
    Rcpp::NumericVector   total_length_bp(N), mean_length(N);

    for (int i = 0; i < N; ++i) {
        group[i]            = fam.samples[i].fid;
        id[i]               = fam.samples[i].iid;
        n_roh[i]            = summaries[i].n_roh;
        n_snps_in_roh[i]    = summaries[i].n_snps_in_roh;
        total_length_bp[i]  = static_cast<double>(summaries[i].total_length_bp);
        mean_length[i]      = (summaries[i].n_roh > 0)
                              ? static_cast<double>(summaries[i].total_length_bp)
                                / summaries[i].n_roh
                              : 0.0;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("group")           = group,
        Rcpp::Named("id")              = id,
        Rcpp::Named("n_ROH")           = n_roh,
        Rcpp::Named("total_length_bp") = total_length_bp,
        Rcpp::Named("mean_length")     = mean_length,
        Rcpp::Named("n_snps_in_roh")   = n_snps_in_roh,
        Rcpp::_["stringsAsFactors"]    = false
    );
}


// ===========================================================================
// C entry point
// ===========================================================================

//' Scan PLINK binary (BED/BIM/FAM) for runs of homozygosity or heterozygosity
//'
//' Low-level C++ entry point called by \code{scanRUNS()}.
//' Do not call directly; use \code{scanRUNS()} instead.
//'
//' @param bed_path Path to the .bed file
//' @param bim_path Path to the .bim file
//' @param fam_path Path to the .fam file
//' @param method Integer: 0 = consecutive (Marras 2015), 1 = sliding window (Bjelland 2013)
//' @param roh_type Integer: 0 = ROHom, 1 = ROHet
//' @param min_snps Minimum SNPs in a qualifying run
//' @param max_opposite Max opposite-type genotypes in a run (consecutive) or window (sliding)
//' @param max_missing Max missing genotypes in a run (consecutive) or window (sliding)
//' @param min_length_bp Minimum run length in base pairs
//' @param max_gap Max gap between consecutive SNPs (bp); >= this breaks a run
//' @param window_size Window width in SNPs (method=1 only)
//' @param threshold Bjelland coverage ratio threshold, strictly > (method=1 only)
//' @param n_threads Number of OpenMP threads (1 = single-threaded)
//' @param verbose If TRUE, print a progress bar during scan and a summary on completion
//' @return Named list: runs, summary, snp_freq, chrom_map
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List C_scan_roh_bed(
    std::string bed_path,
    std::string bim_path,
    std::string fam_path,
    int         method,
    int         roh_type,
    int         min_snps,
    int         max_opposite,
    int         max_missing,
    int         min_length_bp,
    int         max_gap,
    int         window_size,
    double      threshold,
    int         n_threads,
    bool        verbose)
{
    // Parse metadata files
    BimData bim = parse_bim(bim_path);
    FamData fam = parse_fam(fam_path);

    const int n_samples = static_cast<int>(fam.samples.size());
    const int n_snps    = static_cast<int>(bim.snps.size());

    // Build scan parameters
    ScanParams params;
    params.method        = method;
    params.target        = static_cast<int8_t>(roh_type);
    params.min_snps      = min_snps;
    params.max_opposite  = max_opposite;
    params.max_missing   = max_missing;
    params.min_length_bp = min_length_bp;
    params.max_gap       = max_gap;
    params.window_size   = window_size;
    params.threshold     = threshold;
    params.n_threads     = (n_threads > 0) ? n_threads : 1;
    params.verbose       = verbose;

    // Open BED (mmap)
    BedFile bed = open_bed(bed_path, n_samples, n_snps);

    std::vector<RohRecord>    records;
    std::vector<IndivSummary> summaries;
    std::vector<int32_t>      snp_freq;

    try {
        scan_roh_bed(bed, bim, fam, params, records, summaries, snp_freq);
    } catch (...) {
        close_bed(bed);
        throw;
    }
    close_bed(bed);

    // snp_freq as named integer vector
    Rcpp::IntegerVector   r_snp_freq(snp_freq.begin(), snp_freq.end());
    Rcpp::CharacterVector snp_names(n_snps);
    for (int j = 0; j < n_snps; ++j) snp_names[j] = bim.snps[j].snp_name;
    r_snp_freq.attr("names") = snp_names;

    // chrom_map: named integer vector  (chromosome name → uint8_t index)
    const int nc = static_cast<int>(bim.chrom_map.size());
    Rcpp::CharacterVector cmap_names(nc);
    Rcpp::IntegerVector   cmap_idx(nc);
    int k = 0;
    for (const auto& kv : bim.chrom_map) {
        cmap_names[k] = kv.first;
        cmap_idx[k]   = static_cast<int>(kv.second);
        ++k;
    }
    cmap_idx.attr("names") = cmap_names;

    return Rcpp::List::create(
        Rcpp::Named("runs")      = build_runs_df(records, fam, bim),
        Rcpp::Named("summary")   = build_summary_df(summaries, fam),
        Rcpp::Named("snp_freq")  = r_snp_freq,
        Rcpp::Named("chrom_map") = cmap_idx
    );
}


// ===========================================================================
// Binary save / load  (Phase 8)
// ===========================================================================

//' Save ROH scan results to a compact binary file
//'
//' Low-level C++ entry point called by \code{saveROH()}.
//' Do not call directly.
//'
//' @param runs_df   data.frame with columns group, id, chrom, nSNP, from, to, lengthBps
//' @param chrom_map_r Named integer vector from \code{scanRUNS()$chrom_map}
//' @param path      Output file path
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
void C_save_roh(
    Rcpp::DataFrame     runs_df,
    Rcpp::IntegerVector chrom_map_r,
    std::string         path)
{
    // Rebuild chrom_names[256] from the named chrom_map vector
    BimData bim;
    bim.chrom_names.resize(256);
    Rcpp::CharacterVector cn = chrom_map_r.names();
    for (int k = 0; k < static_cast<int>(chrom_map_r.size()); ++k) {
        int idx = chrom_map_r[k];
        if (idx >= 0 && idx < 256)
            bim.chrom_names[idx] = Rcpp::as<std::string>(cn[k]);
    }

    Rcpp::CharacterVector group_col = runs_df["group"];
    Rcpp::CharacterVector id_col    = runs_df["id"];
    Rcpp::CharacterVector chrom_col = runs_df["chrom"];
    Rcpp::IntegerVector   nSNP_col  = runs_df["nSNP"];
    Rcpp::IntegerVector   from_col  = runs_df["from"];
    Rcpp::IntegerVector   to_col    = runs_df["to"];
    const int n = runs_df.nrows();

    // Build unique sample list (first-appearance order) and lookup map
    FamData fam;
    std::unordered_map<std::string, int> sample_idx_map;
    for (int k = 0; k < n; ++k) {
        std::string key = std::string(group_col[k]) + "\t" + std::string(id_col[k]);
        if (sample_idx_map.find(key) == sample_idx_map.end()) {
            sample_idx_map[key] = static_cast<int>(fam.samples.size());
            SampleInfo si;
            si.fid = std::string(group_col[k]);
            si.iid = std::string(id_col[k]);
            fam.samples.push_back(si);
        }
    }

    // Reverse chrom_map: name → uint8 index
    std::unordered_map<std::string, uint8_t> chrom_name_to_idx;
    for (int k = 0; k < static_cast<int>(chrom_map_r.size()); ++k)
        chrom_name_to_idx[Rcpp::as<std::string>(cn[k])] =
            static_cast<uint8_t>(chrom_map_r[k]);

    // Reconstruct RohRecord vector (n_het / n_missing are unavailable from R)
    std::vector<RohRecord> records;
    records.reserve(n);
    for (int k = 0; k < n; ++k) {
        RohRecord rec;
        std::string key = std::string(group_col[k]) + "\t" + std::string(id_col[k]);
        rec.sample_idx = sample_idx_map[key];

        auto it = chrom_name_to_idx.find(Rcpp::as<std::string>(chrom_col[k]));
        rec.chrom_idx  = (it != chrom_name_to_idx.end()) ? it->second : 0;
        rec.start_bp   = from_col[k];
        rec.end_bp     = to_col[k];
        rec.n_snps     = nSNP_col[k];
        rec.n_het      = 0;
        rec.n_missing  = 0;
        records.push_back(rec);
    }

    write_roh_binary(path, records, fam, bim);
}


//' Load ROH scan results from a binary file
//'
//' Low-level C++ entry point called by \code{loadROH()}.
//' Do not call directly.
//'
//' @param path Path to a .roh binary file written by \code{saveROH()}
//' @return Named list with element \code{runs} (data.frame, same format as
//'   \code{scanRUNS()$runs})
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List C_load_roh(std::string path)
{
    RohBinaryData bd = read_roh_binary(path);

    BimData bim;
    bim.chrom_names = bd.chrom_names;

    return Rcpp::List::create(
        Rcpp::Named("runs") = build_runs_df(bd.records, bd.fam, bim)
    );
}


// ===========================================================================
// Permutation-based ROH island detection  (Phase 9)
// ===========================================================================

//' Permutation-based ROH island detection
//'
//' Low-level C++ entry point called by \code{rohIslands()}.
//' Do not call directly; use \code{rohIslands()} instead.
//'
//' @param bed_path    Path to the .bed file
//' @param bim_path    Path to the .bim file
//' @param fam_path    Path to the .fam file
//' @param snp_freq_r  Named integer vector of real SNPROH counts from
//'   \code{scanRUNS()$snp_freq}
//' @param method      Integer: 0=consecutive, 1=sliding
//' @param roh_type    Integer: 0=ROHom, 1=ROHet
//' @param min_snps    Same parameter as the original scan
//' @param max_opposite Same parameter as the original scan
//' @param max_missing Same parameter as the original scan
//' @param min_length_bp Same parameter as the original scan
//' @param max_gap     Same parameter as the original scan
//' @param window_size Sliding window width (sliding method only)
//' @param threshold   Coverage ratio threshold (sliding method only)
//' @param n_threads   OpenMP thread count (parallelism over permutations)
//' @param n_perm      Number of permutations
//' @param percentile  Quantile for threshold derivation (e.g. 0.99)
//' @param seed        MT19937 seed; 0 = draw from random_device
//'
//' @return Named list with \code{thresholds} (named numeric vector, one per
//'   chromosome) and \code{is_island} (named logical vector, one per SNP).
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List C_perm_roh_islands(
    std::string         bed_path,
    std::string         bim_path,
    std::string         fam_path,
    Rcpp::IntegerVector snp_freq_r,
    int                 method,
    int                 roh_type,
    int                 min_snps,
    int                 max_opposite,
    int                 max_missing,
    int                 min_length_bp,
    int                 max_gap,
    int                 window_size,
    double              threshold,
    int                 n_threads,
    int                 n_perm,
    double              percentile,
    int                 seed)
{
    BimData bim = parse_bim(bim_path);
    FamData fam = parse_fam(fam_path);
    const int n_samples = static_cast<int>(fam.samples.size());
    const int n_snps    = static_cast<int>(bim.snps.size());

    ScanParams params;
    params.method        = method;
    params.target        = static_cast<int8_t>(roh_type);
    params.min_snps      = min_snps;
    params.max_opposite  = max_opposite;
    params.max_missing   = max_missing;
    params.min_length_bp = min_length_bp;
    params.max_gap       = max_gap;
    params.window_size   = window_size;
    params.threshold     = threshold;
    params.n_threads     = (n_threads > 0) ? n_threads : 1;
    params.verbose       = false;

    std::vector<int32_t> real_freq(snp_freq_r.begin(), snp_freq_r.end());

    BedFile bed = open_bed(bed_path, n_samples, n_snps);
    PermResult pr;
    try {
        pr = permutation_roh_islands(
            bed, bim, fam, real_freq, params,
            n_perm, percentile,
            static_cast<uint32_t>(seed));
    } catch (...) {
        close_bed(bed);
        throw;
    }
    close_bed(bed);

    // Build named thresholds vector (chromosome name → threshold)
    const int nc = static_cast<int>(pr.chrom_indices.size());
    Rcpp::CharacterVector chrom_names_r(nc);
    Rcpp::NumericVector   thresholds_r(nc);
    for (int k = 0; k < nc; ++k) {
        uint8_t ci = pr.chrom_indices[k];
        chrom_names_r[k] = (static_cast<size_t>(ci) < bim.chrom_names.size() &&
                            !bim.chrom_names[ci].empty())
                           ? bim.chrom_names[ci]
                           : std::to_string(static_cast<int>(ci));
        thresholds_r[k] = pr.thresholds[k];
    }
    thresholds_r.attr("names") = chrom_names_r;

    // Build named is_island logical vector (SNP name → bool)
    Rcpp::LogicalVector   is_island_r(n_snps);
    Rcpp::CharacterVector snp_names_r(n_snps);
    for (int j = 0; j < n_snps; ++j) {
        is_island_r[j] = pr.is_island[j];
        snp_names_r[j] = bim.snps[j].snp_name;
    }
    is_island_r.attr("names") = snp_names_r;

    return Rcpp::List::create(
        Rcpp::Named("thresholds") = thresholds_r,
        Rcpp::Named("is_island")  = is_island_r
    );
}
