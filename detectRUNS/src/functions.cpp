
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
//' @examples
//' data <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA)
//' oppositeAndMissingGenotypes <- findOppositeAndMissing(data, ROHet=TRUE)
//'
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
StringVector findOppositeAndMissing(IntegerVector data, bool ROHet=true) {
  // Initialize oppositeAndMissingGenotypes
  StringVector oppositeAndMissingGenotypes;

  // Initialize vector for names
  std::vector< std::string > names;

  // declare values
  std::string missing = "9";
  std::string opposite = "0";

  // iter in data vector
  for (int i=0; i<data.size(); i++) {
    if (data[i] == NA_INTEGER) {
      // is missing
      oppositeAndMissingGenotypes.push_back(missing);
      // R index are 1 based
      names.push_back(patch::to_string(i+1));
      continue;
    }

    if (ROHet == true){
      if (data[i] == 0) {
        // is homozygote
        oppositeAndMissingGenotypes.push_back(opposite);
        names.push_back(patch::to_string(i+1));
      }
    } else {
      // ROHom condition
      if (data[i] == 1) {
        // is heterozygote
        oppositeAndMissingGenotypes.push_back(opposite);
        names.push_back(patch::to_string(i+1));
      }
    }

  }

  // Finally assign names to StringVector
  oppositeAndMissingGenotypes.attr("names") = names;

  // return results
  return oppositeAndMissingGenotypes;
}


//' Function to slide a window over a vector (individual's genotypes)
//'
//' This is a core function. The functions to detect RUNS are slidden over the genome
//'
//' @param data vector of 0/1/2 genotypes
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
//' @examples #not yet
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
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
//' @param genotype_path genotype (.ped) file location
//'
//' @return a dataframe of POP, ID
//'
//' @examples
//' genotype_path <- system.file("extdata", "subsetChillingham.ped", package = "detectRUNS")
//' pops <- readPOPCpp(genotype_path)
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
//'
// [[Rcpp::export]]
DataFrame readPOPCpp(std::string genotype_path) {
  // the columns of data.frame
  CharacterVector POP;
  CharacterVector ID;

  // open genotype_path
  std::ifstream ifile(genotype_path.c_str());

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
//' @param maxGap max distance between consecutive SNP in a window to be stil considered a potential run
//'
//' @details
//' The consecutive method detect runs by consecutively scanning SNP loci along the genome.
//' No sliding windows are used. Checks on minimum n. of SNP, max n. of opposite and missing genotypes,
//' max gap between adjacent loci and minimum length of the run are implemented (as in the sliding window method).
//' Both runs of homozygosity (RoHom) and of heterozygosity (RoHet) can be search for (option ROHet: TRUE/FALSE)
//'
//' @return A data frame of runs per individual sample
//' @export
//'
//' @examples
//'
// [[Rcpp::export]]
DataFrame consecutiveRunsCpp(IntegerVector indGeno, List individual, DataFrame mapFile,
                             bool ROHet=true, int minSNP=3, int maxOppositeGenotype=1,
                             int maxMiss=1, int minLengthBps=1000, int maxGap=10^6) {

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