#include <Rcpp.h>
#include<string>
#include<map>
using namespace Rcpp;
using namespace std;


//' map two sequences
//' 
//' @param aln gapped sequence
//' @param seq ungapped sequence
//' @return a position map from seq->aln
//' @export
// [[Rcpp::export]]
IntegerVector map_alignment_sequence(std::string & aln, std::string & seq){
  std::map< int, int > m;
  int j=0;
  for(unsigned i=0;i<seq.length();i++){
    if(seq[i]==aln[j]){
      m[i+1]=j+1;
      ++j;
    }else{
      while(aln[j]=='-') ++j;
      m[i+1]=j+1;
      ++j;
    }
  }
  return wrap(m);
}
