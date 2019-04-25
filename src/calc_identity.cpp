#include <Rcpp.h>
#include<string>
#include<map>
#include<vector>
using namespace Rcpp;
using namespace std;

//' calculate identity
//' 
//' @param aln aligned sequences
//' @param len length of aligned sequence
//' @return numeric vector of identical aminoacids between first seq in the aln object and the rest
//' @export 
// [[Rcpp::export]]
NumericVector calc_identity(std::vector<std::string> const &aln, int const len) {
    //int len = aln.size();
    int c=0;
    std::vector<int> mat;
    for(int i=1;i<len;i++){
      c=0;
      for(int j =0;j<aln[i].length();j++){
        if(aln[0][j]==aln[i][j] && aln[0][j]!='-'){
          ++c;
        }
      }
      mat.push_back(c);
    }
    return wrap(mat);
    
}


