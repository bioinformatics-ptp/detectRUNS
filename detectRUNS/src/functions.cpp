#include <Rcpp.h>
#include <string>
#include <sstream>

// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-g
namespace patch {
  template < typename T > std::string to_string( const T& n ) {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

#include <iostream>

using namespace Rcpp;

//' Convert 0/1/2 genotypes to 0/1
//'
//' This is a utility function, that convert 0/1/2 genotypes (AA/AB/BB) into 0/1
//' (either homozygous/heterozygous)
//'
//' @param genotype vector of 0/1/2 genotypes
//'
//' @return converted vector of genotypes (0/1)
//'
//' @examples
//' geno012 <- c(1, 2, 0, 1, NA, 2, 0, NA)
//' geno01 <- genoConvertCpp(geno012)
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
IntegerVector genoConvertCpp(IntegerVector genotype) {
  // deal with missing values (http://adv-r.had.co.nz/Rcpp.html#rcpp-na)
  int missing = NA_INTEGER;

  // deal with map STL class (http://www.cprogramming.com/tutorial/stl/stlmap.html)
  std::map <int, int> converter;

  // set values
  converter[0] = 0;
  converter[1] = 1;
  converter[2] = 0;
  converter[missing] = missing;

  // the converted vector
  IntegerVector results(genotype.size());

  for (int i = 0; i < genotype.size(); i++) {
    results[i] = converter[genotype[i]];
  }

  return results;
}

//' Convert ped genotypes to 0/1
//'
//' This is a utility function, that convert ped genotypes (AA/AB/BB) into 0/1
//' (either homozygous/heterozygous)
//'
//' @param genotype vector of pair of genotypes (01, AA, AG)
//'
//' @return converted vector of genotypes (0/1)
//'
//' @examples
//' ped <- c("A", "A", "A", "B", "5", "5")
//' geno01 <- pedConvertCpp(ped)
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
IntegerVector pedConvertCpp(CharacterVector genotype) {
  // deal with map STL class (http://www.cprogramming.com/tutorial/stl/stlmap.html)
  std::map <std::string, int> missing;
  std::string allele1, allele2;

  // Loading message function from R (https://github.com/RcppCore/Rcpp/issues/195)
  Rcpp::Function warn("warning");
  Rcpp::Function stop("stop");

  // check genotype size
  if (genotype.size() % 2 != 0 ) {
    stop(std::string("Need .ped input with 2 alleles per marker"));
  }

  // set values
  missing["0"] = NA_INTEGER;
  missing["5"] = NA_INTEGER;
  missing["N"] = NA_INTEGER;

  // the converted vector
  IntegerVector results(genotype.size() / 2);

  for (int i = 0; i < genotype.size(); i=i+2) {
    // for semplicity
    allele1 = genotype[i];
    allele2 = genotype[i+1];

    if (allele1 != allele2) {
      // test for only one allele missing (http://www.cprogramming.com/tutorial/stl/stlmap.html)
      // check that one allele isn't in the missing map structure
      if (missing.find(allele1) != missing.end() || missing.find(allele2) != missing.end()) {
        warn(std::string("Found only one allele missing in a pair"));
      }

      // Set allele as heterozygous
      results[i/2] = 1;

    } else {
      //test for missing allele (allele are equals condition)
      if (missing.find(allele1) != missing.end()) {
        results[i/2] = NA_INTEGER;

      } else {
        // alleles are homozygous
        results[i/2] = 0;
      }

    } //allele are equals condition

  } // cicle for allele pair

  return results;
}



//' Function to check whether a window is (loosely) homozygous or not
//'
//' This is a core function. Parameters on how to consider a window homozygous are here (maxHet, maxMiss)
//'
//' @param x vector of 0/1 genotypes (from genoConvert())
//' @param gaps vector of differences between consecutive positions (gaps) in bps
//' @param maxHet max n. of heterozygous SNP in a homozygous window
//' @param maxMiss max n. of missing in a window
//' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
//'
//' @return TRUE/FALSE (whether a window is homozygous or NOT)
//'
//' @examples
//' maxHom <- 1
//' maxMiss <- 1
//' maxGap <- 10^6
//' x <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
//'        0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' gaps <- c(3721, 3871, 7059, 4486, 7545, 4796, 3043, 9736, 3495, 5051,
//'           9607, 6555, 11934, 6410, 3415, 1302, 3110, 6609, 3292)
//' test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
//' # test is true
//' x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
//'        1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
//'           615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)
//' test <- homoZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
//' # test is false
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
bool homoZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHet, int maxMiss, int maxGap) {
  // check gaps
  for (int i=0; i< gaps.size(); i++) {
    if (gaps[i] > maxGap) {
      return false;
    }
  }

  // count Heterozygots
  int nHet = std::count(x.begin(), x.end(), 1);

  // count missing values
  int nMiss = std::count(x.begin(), x.end(), NA_INTEGER);

  if (nHet > maxHet || nMiss > maxMiss) {
    return false;
  }

  // if I pass all checks
  return true;
}


//' Function to check whether a window is (loosely) heterozygous or not
//'
//' This is a core function. Parameters on how to consider a window heterozygous are here (maxHom, maxMiss)
//'
//' @param x vector of 0/1 genotypes (from genoConvert())
//' @param gaps vector of differences between consecutive positions (gaps) in bps
//' @param maxHom max n. of homozygous SNP in a heterozygous window
//' @param maxMiss max n. of missing in a window
//' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
//'
//' @return TRUE/FALSE (whether a window is heterozygous or NOT)
//' @export
//'
//' @examples
//' maxHom <- 1
//' maxMiss <- 1
//' maxGap <- 10^6
//' x <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
//'        1, 1, 1, 1, 0, 0, 1, 0, 0, 0)
//' gaps <- c(4374, 8744, 5123, 14229, 5344, 690, 8566, 5853, 2369, 3638,
//'           4848, 600, 2333, 976, 2466, 2269, 5411, 6021, 4367)
//' test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
//' # test is false
//' x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
//'        1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' gaps <- c(2514, 2408, 2776, 2936, 1657, 494, 1436, 680, 909, 678,
//'           615, 1619, 2058, 2446, 1085, 660, 1259, 1042, 2135)
//' test <- heteroZygotTestCpp(x, gaps, maxHom, maxMiss, maxGap)
//' # test is true
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
bool heteroZygotTestCpp(IntegerVector x, IntegerVector gaps, int maxHom, int maxMiss, int maxGap) {
  // check gaps
  for (int i=0; i< gaps.size(); i++) {
    if (gaps[i] > maxGap) {
      return false;
    }
  }

  // count Homozygots
  int nHom = std::count(x.begin(), x.end(), 0);

  // count missing values
  int nMiss = std::count(x.begin(), x.end(), NA_INTEGER);

  if (nHom > maxHom || nMiss > maxMiss) {
    return false;
  }

  // if I pass all checks
  return true;
}


//' Function to slide a window over a vector (individual's genotypes)
//'
//' This is a core function. The functions to detect RUNS are slidden over the genome
//'
//' @param data vector of pair of genotypes (01, AA, AG)
//' @param gaps vector of differences between consecutive positions (gaps) in bps
//' @param windowSize size of window (n. of SNP)
//' @param step by which (how many SNP) is the window slidden
//' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
//' @param ROHet shall we detect ROHet or ROHom?
//' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
//' @param maxMiss max. n. of missing SNP
//'
//' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
//'
//' @examples #not yet
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
LogicalVector slidingWindowCpp(CharacterVector data, IntegerVector gaps, int windowSize, int step,
                               int maxGap, bool ROHet=true, int maxOppositeGenotype=1, int maxMiss=1) {

  // Loading message function from R (https://github.com/RcppCore/Rcpp/issues/195)
  Rcpp::Function msg("message");

  // convert genotype
  IntegerVector y = pedConvertCpp(data);

  // get data lenght
  int data_length = y.size();

  // calculate spots size
  int spots_lenght = (data_length - windowSize) / step +1;

  // initialize results
  LogicalVector results(spots_lenght, false);

  // declare iterators
  IntegerVector::const_iterator from, to;

  // eval RoHet or RoHom
  // if (ROHet == true) {
  //   msg(std::string("Analysing Runs of Heterozygosity (ROHet)"));
  // } else {
  //   msg(std::string("Analysing Runs of Homozygosity (ROHom)"));
  // }

  // evaluating windows
  for (int i=0; i<spots_lenght; i++) {
    // calculate y_spots
    from = y.begin() + i*step;
    to = y.begin() + i*step + windowSize;

    //Slice the original vector
    IntegerVector y_spots(from, to);

    // calculate gaps_spots
    from = gaps.begin() + i*step;
    to = gaps.begin() + i*step + windowSize -1;
    IntegerVector gaps_spots(from, to);

    // debug
    // Rcout << y_spots << std::endl;
    // Rcout << gaps_spots << std::endl;

    // eval RoHet or RoHom
    if (ROHet == true) {
      // calculate result
      results[i] = heteroZygotTestCpp(y_spots, gaps_spots, maxOppositeGenotype, maxMiss, maxGap);
      // Rcout << results[i] << std::endl;

    } else {
      // calculate result
      results[i] = homoZygotTestCpp(y_spots, gaps_spots, maxOppositeGenotype, maxMiss, maxGap);
      // Rcout << results[i] << std::endl;
    }

  }

  // check this affermation (could be N of SNPs - window +1)
  // msg(std::string("Length of homozygous windows overlapping SNP loci (should be equal to the n. of SNP in the file): "), results.size());

  return results;
}


//' Function to return a vector of T/F for whether a SNP is or not in a RUN
//'
//' This is a core function. The function to determine whether a SNP is or not in a RUN.
//' The ratio between homozygous/heterozygous windows and total n. of windows is computed here
//'
//' @param RunVector vector of TRUE/FALSE (is a window homozygous/heterozygous?)
//' @param windowSize size of window (n. of SNP)
//' @param threshold threshold to call a SNP in a RUN
//'
//' @return vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
//'
//' @examples #not yet
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
LogicalVector snpInRunCpp(LogicalVector RunVector, const int windowSize, const float threshold) {
  // Loading message function from R (https://github.com/RcppCore/Rcpp/issues/195)
  Rcpp::Function msg("message");

  // get vector size
  int RunVector_length = RunVector.size();

  // debug
  // msg(std::string("Length of input vector: "), RunVector_length);
  // msg(std::string("Window size: "), windowSize);
  // msg(std::string("Threshold for calling SNP in a Run: ") + patch::to_string(threshold));

  // compute total n. of overlapping windows at each SNP locus (see Bjelland et al. 2013)
  // initialize a vector with window as default value
  std::vector<int> nWin(RunVector_length, windowSize);

  // then fix values at both sides
  for (int i=0; i<windowSize; i++) {
    nWin[i] = i+1;
    nWin[RunVector_length - i-1] = i+1;
  }

  // compute n. of homozygous/heterozygous windows that overlap at each SNP locus (Bjelland et al. 2013)
  float hWin;
  float quotient;
  LogicalVector::iterator from, to;

  // the returned value, a logical vector with false values as default values
  LogicalVector snpRun(RunVector_length, false);

  for (int i=0; i < RunVector_length; i++) {
    //get from-to index fom nWin. Get iterators
    from = RunVector.begin() + i;
    // the to index iterator is excluded
    to = RunVector.begin() + i + nWin[i];

    //count TRUE in interval
    hWin = std::count(from, to, true);

    //calc quotient
    quotient = hWin/nWin[i];

    //vector of SNP belonging to a ROH. True if yes (quotient > threshold)
    if (quotient > threshold) {
      snpRun[i] = true;
    }

    //debug
    // if (i > 20 && i < 30) {
    //   Rcout << "i: " << i << " hWin: " << hWin << " nWin[i]: "<< nWin[i] << " quotient: " << quotient;
    //   Rcout << " snpRun[i]: " << snpRun[i] << std::endl;
    // }

  }

  //debug
  //msg(std::string("Lenght of output vector: "), snpRun.size());

  return snpRun;
}
