
#include <Rcpp.h>
#include <string>
#include <sstream>
#include <fstream>


// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-g
namespace patch {
  template < typename T > std::string to_string( const T& n ) {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


// http://stackoverflow.com/questions/17694579/use-stdfill-to-populate-vector-with-increasing-numbers
struct IncGenerator {
  int current_;
  IncGenerator (int start) : current_(start) {}
  int operator() () { return current_++; }
};


#include <iostream>


using namespace Rcpp;


// http://gallery.rcpp.org/articles/fast-factor-generation/
template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "factor";
  return out;
}

// [[Rcpp::export]]
SEXP fast_factor( SEXP x ) {
  switch( TYPEOF(x) ) {
  case INTSXP: return fast_factor_template<INTSXP>(x);
  case REALSXP: return fast_factor_template<REALSXP>(x);
  case STRSXP: return fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}


//' Convert 0/1/2 genotypes to 0/1
//'
//' This is a utility function, that convert 0/1/2 genotypes (AA/AB/BB) into 0/1
//' (either homozygous/heterozygous)
//'
//' @param genotype vector of 0/1/2 genotypes
//'
//' @return converted vector of genotypes (0/1)
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
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
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
IntegerVector pedConvertCpp(CharacterVector genotype) {
  // deal with map STL class (http://www.cprogramming.com/tutorial/stl/stlmap.html)
  std::map <std::string, int> missing;
  std::string allele1, allele2;

  // check genotype size
  if (genotype.size() % 2 != 0 ) {
    stop(std::string("Need .ped input with 2 alleles per marker"));
  }

  // set all possible missing values
  missing["0"] = NA_INTEGER;
  missing["5"] = NA_INTEGER;
  missing["N"] = NA_INTEGER;
  missing["-"] = NA_INTEGER;

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
        stop(std::string("Found only one allele missing in a pair"));
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
//' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
//'
//' @return TRUE/FALSE (whether a window is homozygous or NOT)
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
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
//' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
//'
//' @return TRUE/FALSE (whether a window is heterozygous or NOT)
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
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


//' Function to calculate oppositeAndMissingGenotypes array
//'
//' This is an helper function, this will be called by another function
//'
//' @param data vector of 0/1/2 genotypes
//' @param ROHet TRUE in ROHet evaluation, FALSE for ROHom
//'
//' @return character array; names will be index in which opposite and missing
//' snps are found in data array
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
StringVector findOppositeAndMissing(IntegerVector data, bool ROHet=true) {
  // Initialize a temporary array for oppositeAndMissingGenotypes
  std::vector<std::string> tmp(data.size());

  // Initialize vector for names
  std::vector<std::string> names(data.size());

  // declare values
  std::string missing = "9";
  std::string opposite = "0";

  // count missing and opposite
  int index = 0;

  // iter in data vector
  for (int i=0; i<data.size(); i++) {
    if (data[i] == NA_INTEGER) {
      // is missing
      tmp[index]= missing;
      // R index are 1 based
      names[index] = patch::to_string(i+1);
      // increment index
      index++;
      continue;
    }

    if (ROHet == true){
      if (data[i] == 0) {
        // is homozygote
        tmp[index] = opposite;
        names[index] = patch::to_string(i+1);
        // increment index
        index++;
      }
    } else {
      // ROHom condition
      if (data[i] == 1) {
        // is heterozygote
        tmp[index] = opposite;
        names[index] = patch::to_string(i+1);
        // increment index
        index++;
      }
    }

  }

  // resize StringVectors
  tmp.resize(index);
  names.resize(index);

  // get a string vector from tmp data
  StringVector oppositeAndMissingGenotypes = wrap(tmp);

  // Finally assign names to StringVector
  oppositeAndMissingGenotypes.attr("names") = names;

  // return results
  return oppositeAndMissingGenotypes;
}


//' Function to slide a window over a vector (individual's genotypes)
//'
//' This is a core function. The functions to detect RUNS are slid over the genome
//'
//' @param data vector of 0/1/2 genotypes
//' @param gaps vector of differences between consecutive positions (gaps) in bps
//' @param windowSize size of window (n. of SNP)
//' @param step by which (how many SNP) is the window slid
//' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
//' @param ROHet shall we detect ROHet or ROHom?
//' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
//' @param maxMiss max. n. of missing SNP
//'
//' @return vector of TRUE/FALSE (whether a window is homozygous or NOT)
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
List slidingWindowCpp(IntegerVector data, IntegerVector gaps, int windowSize,
                      int step, int maxGap, bool ROHet=true,
                      int maxOppositeGenotype=1, int maxMiss=1) {

  // get data lenght
  int data_length = data.size();

  // calculate spots size
  int spots_lenght = (data_length - windowSize) / step +1;

  // initialize results
  LogicalVector windowStatus(spots_lenght, false);

  // calculate opposite and missing snps
  StringVector oppositeAndMissingGenotypes = findOppositeAndMissing(data, ROHet);

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
    from = data.begin() + i*step;
    to = data.begin() + i*step + windowSize;

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
      windowStatus[i] = heteroZygotTestCpp(y_spots, gaps_spots, maxOppositeGenotype, maxMiss, maxGap);
      // Rcout << windowStatus[i] << std::endl;

    } else {
      // calculate result
      windowStatus[i] = homoZygotTestCpp(y_spots, gaps_spots, maxOppositeGenotype, maxMiss, maxGap);
      // Rcout << windowStatus[i] << std::endl;
    }

  }

  // check this affermation (could be N of SNPs - window +1)
  // msg(std::string("Length of homozygous windows overlapping SNP loci (should be equal to the n. of SNP in the file): "), windowStatus.size());

  // define a list of results and return it
  return List::create(Named("windowStatus")=windowStatus,
                      Named("oppositeAndMissingGenotypes")=oppositeAndMissingGenotypes);
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
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
LogicalVector snpInRunCpp(LogicalVector RunVector, const int windowSize, const float threshold) {
  // get vector size
  int RunVector_length = RunVector.size();

  // compute total n. of overlapping windows at each SNP locus (see Bjelland et al. 2013)
  // initialize a vector with window as default value. I need to have nWin size as
  // the number of windows is size would be 1
  int nWin_length = RunVector_length + windowSize -1 ;
  std::vector<int> nWin(nWin_length, windowSize);

  // then fix values at both sides
  for (int i=0; i<windowSize; i++) {
    nWin[i] = i+1;
    nWin[nWin_length - i-1] = i+1;
  }

  // create two sets of indices to slice the vector of windows containing or not
  // a run (RunVector). Since they are array of positions, they should refere to
  // O based coordinates, since in C++ array starts at 0 index
  std::vector<int> ind1(nWin_length, 0);

  // initialize with the first element of array.
  IncGenerator g1(0);

  // Fill with the result of calling g1() repeatedly.
  // ind1 = c(rep(1,windowSize-1),seq(1,RunVector_length))
  // since ind1 is a vector of indexes, need to numbering after than R
  std::generate(ind1.begin()+windowSize-1, ind1.end(), g1);

  // define ind2 vector
  std::vector<int> ind2(nWin_length, RunVector_length-1);

  // initialize with the first element of array
  IncGenerator g2(0);

  // Fill with the result of calling g2() repeatedly.
  // ind2 = c(seq(1,RunVector_length),rep(RunVector_length,windowSize-1))
  std::generate(ind2.begin(), ind2.end()-windowSize+1, g2);

  // compute n. of homozygous/heterozygous windows that overlap at each SNP locus (Bjelland et al. 2013)
  float hWin;
  float quotient;
  LogicalVector::iterator from, to;

  // the returned value, a logical vector with false values as default values
  LogicalVector snpRun(nWin_length, false);

  for (int i=0; i < nWin_length; i++) {
    // get from-to index fom nWin. Get iterators
    from = RunVector.begin() + ind1[i];

    // the to index iterator is excluded
    to = RunVector.begin() + ind2[i] + 1;

    // count TRUE in interval
    hWin = std::count(from, to, true);

    //calc quotient
    quotient = hWin/nWin[i];

    //vector of SNP belonging to a ROH. True if yes (quotient > threshold)
    if (quotient > threshold) {
      snpRun[i] = true;
    }

  }

  return snpRun;
}


//' Function to return a dataframe of population (POP, ID)
//'
//' This is a core function. Read PED file and returns a data.frame with the first two
//' columns
//'
//' @param genotypeFile genotype (.ped) file location
//'
//' @return a dataframe of POP, ID
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
DataFrame readPOPCpp(std::string genotypeFile) {
  // the columns of data.frame
  CharacterVector POP;
  CharacterVector ID;

  // open genotypeFile
  std::ifstream ifile(genotypeFile.c_str());

  // we read the full line here
  std::string line;

  // A tocken for read columns
  std::string token;

  // read the current line (http://stackoverflow.com/questions/30181600/reading-two-columns-in-csv-file-in-c)
  while (std::getline(ifile, line)) {
    // construct a string stream from line
    std::istringstream iss(line);

    // current token
    std::getline(iss, token, ' ');
    POP.push_back(token);

    std::getline(iss, token, ' ');
    ID.push_back(token);
  }

  // This is the resulting new dataframe (New Data Frame). A Cpp instance of
  // stringsAsFactors = TRUE data.frame
  DataFrame NDF = DataFrame::create(Named("POP")=POP, Named("ID")=ID, _["stringsAsFactors"] = false);

  return NDF;
}


// helper RunData structure for consecutiveRunsCpp
struct RunData {
  int nOpposite;
  int nMiss;
  int runH;
  int lengte;
  std::string chrom;
  int start;
  int end;
};


// helper initializeRun function for consecutiveRunsCpp
RunData initializeRun(std::string chrom, int start) {
  RunData run_data;

  // initialize values
  run_data.chrom = chrom;
  run_data.start = start;
  run_data.end = start;
  run_data.nOpposite = 0;
  run_data.nMiss = 0;
  run_data.runH = 0;
  run_data.lengte = 0;

  return run_data;
}

// helper function to update RUNs data for consecutiveRunsCpp
void updateRUNs(RunData run_data, std::string iid, std::string fid, CharacterVector *group,
                CharacterVector *id, CharacterVector *chrom, IntegerVector *nSNP,
                IntegerVector *from, IntegerVector *to, IntegerVector *lengthBps) {

  // update arrays for current RUN
  group->push_back(fid);
  id->push_back(iid);
  chrom->push_back(run_data.chrom);
  nSNP->push_back(run_data.runH);
  from->push_back(run_data.start);
  to->push_back(run_data.end);
  lengthBps->push_back(run_data.lengte);

}

//' Function to detect consecutive runs in a vector (individual's genotypes)
//'
//' This is a core function. It implements the consecutive method for detection of runs in diploid genomes
//' (see Marras et al. 2015)
//'
//' @param indGeno vector of 0/1/NAs of individual genotypes (0: homozygote; 1: heterozygote)
//' @param individual list of group (breed, population, case/control etc.) and ID of individual sample
//' @param mapFile Plink map file (for SNP position)
//' @param ROHet shall we detect ROHet or ROHom?
//' @param minSNP minimum number of SNP in a run
//' @param maxOppositeGenotype max n. of homozygous/heterozygous SNP
//' @param maxMiss max. n. of missing SNP
//' @param minLengthBps min length of a run in bps
//' @param maxGap max distance between consecutive SNP in a window to be still considered a potential run
//'
//' @details
//' The consecutive method detect runs by consecutively scanning SNP loci along the genome.
//' No sliding windows are used. Checks on minimum n. of SNP, max n. of opposite and missing genotypes,
//' max gap between adjacent loci and minimum length of the run are implemented (as in the sliding window method).
//' Both runs of homozygosity (RoHom) and of heterozygosity (RoHet) can be search for (option ROHet: TRUE/FALSE)
//'
//' @return A data frame of runs per individual sample
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
DataFrame consecutiveRunsCpp(IntegerVector indGeno, List individual, DataFrame mapFile,
                             bool ROHet=true, int minSNP=3, int maxOppositeGenotype=1,
                             int maxMiss=1, int minLengthBps=1000, int maxGap=10e5) {

  // Bool to int conversion: ifelse(ROHet,1,0)
  int typ = ROHet;

  // animal data lile IID or FID
  std::string iid = individual["IID"];
  std::string fid = individual["FID"];

  // initialize variables. First chromosome in ordered mapFile
  CharacterVector Chrom = mapFile["Chrom"];
  IntegerVector bps = mapFile["bps"];

  std::string lastChrom = as<std::string>(Chrom[0]);
  int lastPos = bps[0];
  std::string currentChrom;
  int currentPos;
  int gap;

  // a flag to determine if I'm in a RUN or not
  bool flag_run = false;
  RunData run_data;

  // the columns of data.frame Defining data types accordingly slinding window
  CharacterVector group;
  CharacterVector id;
  CharacterVector chrom;
  IntegerVector nSNP;
  IntegerVector from;
  IntegerVector to;
  IntegerVector lengthBps;

  for (int i = 0; i < indGeno.size(); i++) {
    // Check for Chromosome
    currentChrom = Chrom[i];

    // get current position
    currentPos = bps[i];

    // test if chromosome is changed
    if (currentChrom != lastChrom ) {
      // if I have run, write to file
      if(flag_run == true && run_data.runH >= minSNP && run_data.lengte >= minLengthBps) {
        // debug
        // Rcout << "Update RUN: chromosome changed" << std::endl;
        updateRUNs(run_data, iid, fid, &group, &id, &chrom, &nSNP, &from, &to, &lengthBps);
      }

      // update chrom and positions
      lastChrom = currentChrom;
      lastPos = currentPos;

      // unset RUN flag. New runs with new chromosomes!!!
      flag_run = false;
    }

    // calculate gap between consecutive SNP (in the same chrom)
    gap = currentPos - lastPos;

    // check if current gap is larger than max allowed gap. No matter current SNP
    if (gap >= maxGap) {
      if(flag_run == true && run_data.runH >= minSNP && run_data.lengte >= minLengthBps) {
        // debug
        // Rcout << "Update RUN: gap size exceeded" << std::endl;
        updateRUNs(run_data, iid, fid, &group, &id, &chrom, &nSNP, &from, &to, &lengthBps);
      }

      // unset RUN flag
      flag_run = false;
      }

    // Start for ==. Is a new RUN or not? ensure that indGeno[i] is a number
    if (indGeno[i] == typ && indGeno[i] != NA_INTEGER){
      // initialize run if not yet initialized, or just written after a big GAP
      if (flag_run == false) {
        //debug
        // Rcout << "Creating new RUN at i = " << i << std::endl;
        // run_data is a struct of attributes for the current RUN
        run_data = initializeRun(currentChrom, currentPos);
        flag_run = true;
      }

      // update run_data values
      run_data.runH++;
      run_data.end = currentPos;
      run_data.lengte = (run_data.end - run_data.start);

    } // condition: the genotype I want

    // start if !=
    else if (indGeno[i] != typ && indGeno[i] != NA_INTEGER){
      // if not in a run, don't do anything
      if (flag_run == false) {
        continue;
      }

      // update nOpposite genotypes (in a run)
      run_data.nOpposite++;

      // check if maxOppositeGenotype is reached
      if (run_data.nOpposite <= maxOppositeGenotype) {
      // update run_data values. This opposite genotype is a part of the RUN
        run_data.runH++;
        run_data.end = currentPos;
        run_data.lengte = (run_data.end - run_data.start);

      } else {
        // debug
        // Rcout << "max opposite reached" << std::endl;

        if (run_data.runH >= minSNP && run_data.lengte >= minLengthBps) {
          updateRUNs(run_data, iid, fid, &group, &id, &chrom, &nSNP, &from, &to, &lengthBps);
        }

        // unset RUN flag
        flag_run = false;

      } // condition nOpposite greather than maxOppositeGenotype

    } // condition: opposite genotype

    // start if 'NA'
    else if (indGeno[i] == NA_INTEGER) {
      // if not in a run, don't do anything
      if (flag_run == false) {
        continue;
      }

      // update missing values
      run_data.nMiss++;

      // check if maxMiss is reached
      if (run_data.nMiss <= maxMiss){
        // update run_data values. This missing genotype is a part of the RUN
        run_data.runH++;
        run_data.end = currentPos;
        run_data.lengte = (run_data.end - run_data.start);
      } else {
        // debug
        // Rcout << "max missing reached" << std::endl;

        if (run_data.runH >= minSNP && run_data.lengte >= minLengthBps) {
          updateRUNs(run_data, iid, fid, &group, &id, &chrom, &nSNP, &from, &to, &lengthBps);
        }

        // unset RUN flag
        flag_run = false;

      } // condition: nMissing greather than permitted

    } // condition: missing genotype

  // update positions
  lastPos = currentPos;

  } // cicle: for i in genotypes

  // last snp if it is in a run
  if (flag_run == true ) {
    if (run_data.runH >= minSNP && run_data.lengte >= minLengthBps) {
      // Rcout << "Last RUN finished with last SNP" << std::endl;
      updateRUNs(run_data, iid, fid, &group, &id, &chrom, &nSNP, &from, &to, &lengthBps);
    }

    // unset RUN flag
    flag_run = false;
  }

  // initialize dataframe of results.
  DataFrame res = DataFrame::create(
    Named("group")=group, Named("id")=id, Named("chrom")=chrom, Named("nSNP")=nSNP,
    Named("from")=from, Named("to")=to, Named("lengthBps")=lengthBps,
    _["stringsAsFactors"] = false);

  // debug
  if(res.nrows() > 0) {
    Rcout << "N. of RUNS for individual " << iid << " is: " << res.nrows() << std::endl;
  } else {
    Rcout << "No RUNs found for animal " << iid << std::endl;
  }

  // returning all runs for this individual genotype
  return(res);

}


// Helper class to deal with runs
class Runs {
  std::vector<std::string> breed;
  std::vector<std::string> chromosome;
  std::vector<int> start;
  std::vector<int> end;
  int size;

  // a map object for position indexing: in order to avoid to scan all runs
  // for a SNP, I will divide all chromosome size in chunk and in each chunk
  // I will put the run index which span this chunk. When I will have a position
  // I will calculate its chunk and then I will scan for each runs in this chunk
  // if the snps belong to it or not
  std::map <int, std::vector<int> > indexes;

  // chunk by size
  int chunk;

public:
  Runs(DataFrame runs);
  // function to count runs by breed
  std::map <std::string, int> countSnpbyBreed(
      int position, std::vector<std::string> unique_breeds);
  void dumpRuns();
};

Runs::Runs(DataFrame runs) {
  // declare variables
  int start, end;

  // define chunk size
  this->chunk = 1e6;

  // get vectors for simplicity
  this->breed = as<std::vector<std::string> >(runs["POPULATION"]);
  this->chromosome = as<std::vector<std::string> >(runs["CHROMOSOME"]);
  this->start = as<std::vector<int> >(runs["START"]);
  this->end = as<std::vector<int> >(runs["END"]);

  // set size
  this->size = this->chromosome.size();

  // index the object. Get index position by dividing run position by chunks
  for (int i=0; i<this->size; i++) {
    // by dividing integer number, I get the integer part of division
    start = this->start[i] / this->chunk;
    end = this->end[i] / this->chunk;

    // add this run index to indexes map
    for (int j=start; j<=end; j++){
      if (indexes.find(j) == indexes.end()) {
        indexes[j] = std::vector<int>();
      }

      // push back index in chunk list
      indexes[j].push_back(i);
    }
  }
}

// count how many type a snp (position) belong to a RUN
std::map <std::string, int> Runs::countSnpbyBreed(
    int position, std::vector<std::string> unique_breeds) {
  // define variables
  std::map <std::string, int> counts;
  std::string ras;

  // initialize
  for (unsigned int i=0; i<unique_breeds.size(); i++) {
    ras = unique_breeds[i];
    counts[ras] = 0;
  }

  // calculating chunk position: we will investigate only runs spanning this region
  int index = position / this->chunk;

  // debug
  // Rcout << "Got position: " << position << " (";
  // Rcout << index * this->chunk << "-" << (index+1) * this->chunk;
  // Rcout << ")" << std::endl;

  // get chunk indexes
  std::vector<int> chunk_indexes = this->indexes[index];

  // iter over runs in chunks
  for(std::vector<int>::iterator it = chunk_indexes.begin();
      it != chunk_indexes.end(); ++it) {

    if (position >= this->start[*it] && position <= this->end[*it]) {
      counts[this->breed[*it]]++;
    }
  }

  return counts;
}

void Runs::dumpRuns() {
  for (int i=0; i<this->size; i++) {
    Rcout << "breed " << this->breed[i] << " chrom " << this->chromosome[i];
    Rcout << " start " << this->start[i] << " end " << this->end[i] << std::endl;
  }

  // dump indexes. Cicle around map iterators
  std::map <int, std::vector<int> >::iterator it;

  for (it = this->indexes.begin(); it != this->indexes.end(); ++it) {
    Rcout << "Chunk " << it->first * this->chunk << "-" << (it->first+1) * this->chunk << ": ";
    for(std::vector<int>::iterator it2 = it->second.begin();
        it2 != it->second.end(); ++it2) {
      Rcout << *it2+1 << " ";
    }

    Rcout << std::endl;
  }
}


//' Function to count number of times a SNP is in a RUN
//'
//'
//' @param runsChrom R object (dataframe) with results per chromosome
//' @param mapChrom R map object with SNP per chromosome
//' @param genotypeFile genotype (.ped) file location
//'
//' @return dataframe with counts per SNP in runs (per population)
//'
//' @import utils
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
DataFrame snpInsideRunsCpp(DataFrame runsChrom, DataFrame mapChrom,
                           std::string genotypeFile) {

  // transform R object in Cpp object
  std::vector<int> POSITIONS = as<std::vector<int> >(mapChrom["POSITION"]);
  std::vector<std::string> SNP_NAME = as<std::vector<std::string> >(mapChrom["SNP_NAME"]);
  std::vector<int> CHR = as<std::vector<int> >(mapChrom["CHR"]);

  // get unique breeds
  CharacterVector population = runsChrom["POPULATION"];
  std::vector<std::string> unique_breeds = as<std::vector<std::string> >(unique(population));

  // sort unique breeds
  std::sort(unique_breeds.begin(), unique_breeds.end());

  // declare others variables
  std::string ras;
  int pos, index, map_size = SNP_NAME.size();
  std::map <std::string, int> snpCounts, nBreeds ;

  // define result size like n SNPs * unique_breeds
  int result_size = map_size * unique_breeds.size();

  // the columns of data.frame Defining data types accordingly slinding window
  CharacterVector snp_name(result_size);
  IntegerVector chrom(result_size); // as defined in R function. What about X?
  IntegerVector position(result_size);
  IntegerVector count(result_size);
  CharacterVector breed(result_size);
  NumericVector percentage(result_size);

  // get all populations
  DataFrame pops = readPOPCpp(genotypeFile);
  CharacterVector pop = pops["POP"];

  // instantiate a Runs object
  Runs runs(runsChrom);

  // debug: dump object
  // runs.dumpRuns();

  // cicle among single breeds and find breed numbers
  for (unsigned int i=0; i<unique_breeds.size(); i++) {
    ras = unique_breeds[i];

    // get total if individuals by breed
    nBreeds[ras] = std::count(pop.begin(), pop.end(), ras.c_str());
    // Rcout << "N. of animals of Population " << ras << ": " << nBreeds[ras] << std::endl;
  }

  // iterate over position. Update single values
  for (unsigned int j=0; j<POSITIONS.size(); j++) {
    // get a snp position
    pos = POSITIONS[j];

    // update counts
    snpCounts = runs.countSnpbyBreed(pos, unique_breeds);

    // update results by breed
    for (unsigned int i=0; i<unique_breeds.size(); i++) {
      // get a breed
      ras = unique_breeds[i];

      //calculating results index
      index = i * map_size + j;

      // update values
      count[index] = snpCounts[ras];
      breed[index] = ras;
      percentage[index] = double(snpCounts[ras])/nBreeds[ras]*100;

      // debug
      // if (SNP_NAME[j] == "OAR24_6970428.1") {
      //   Rcout << "Index i: " << i << " Index j: " << j << " Final index: " << index;
      //   Rcout << " Breed: " << ras << " Count: " << snpCounts[ras] << std::endl;
      // }

      // read data from mapChrom
      snp_name[index] = SNP_NAME[j];
      chrom[index] = CHR[j];
      position[index] = pos;

    } // cicle for snp position

  } // cicle for breed

  // initialize dataframe of results.
  DataFrame res = DataFrame::create(
    Named("SNP_NAME")=snp_name, Named("CHR")=chrom, Named("POSITION")=position,
    Named("COUNT")=count, Named("BREED")=fast_factor(breed)  , //returns a factor
    Named("PERCENTAGE")=percentage, _["stringsAsFactors"] = false);

  // returning all runs for this individual genotype
  return(res);
}
