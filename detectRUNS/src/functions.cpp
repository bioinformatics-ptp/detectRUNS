#include <Rcpp.h>
using namespace Rcpp;

//' Function to return a vector of T/F for whether a SNP is or not in a RUN
//'
//' This is a core function. The function to determine whether a SNP is or not in a RUN.
//' The ratio between homozygous/heterozygous windows and total n. of windows is computed here
//'
//' @param RunVector vector of TRUE/FALSE (is a window homozygous/heterozygous?)
//' @param window size of window (n. of SNP)
//' @param threshold threshold to call a SNP in a RUN
//'
//' @return vector of TRUE/FALSE (whether a SNP is in a RUN or NOT)
//'
//' @examples #not yet
//' @useDynLib detectRUNS
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector snpInRunCpp(LogicalVector RunVector, const int window, const float threshold) {
  // get vector size
  int RunVector_length = RunVector.size();

  Rcout << "Length of input vector: " << RunVector_length << std::endl;
  Rcout << "Window size: " << window << std::endl;
  Rcout << "Threshold for calling SNP in a Run: " << threshold << std::endl;

  // compute total n. of overlapping windows at each SNP locus (see Bjelland et al. 2013)
  // initialize a vector with window as default value
  std::vector<int> nWin(RunVector_length, window);

  // then fix values at both sides
  for (int i=0; i<window; i++) {
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
  Rcout << "Lenght of output file: " << snpRun.size() << std::endl;

  return snpRun;
}
