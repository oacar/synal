// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calc_identity
NumericVector calc_identity(std::vector<std::string> const& aln, int const len);
RcppExport SEXP _synal_calc_identity(SEXP alnSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type aln(alnSEXP);
    Rcpp::traits::input_parameter< int const >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_identity(aln, len));
    return rcpp_result_gen;
END_RCPP
}
// map_alignment_sequence
IntegerVector map_alignment_sequence(std::string& aln, std::string& seq);
RcppExport SEXP _synal_map_alignment_sequence(SEXP alnSEXP, SEXP seqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type aln(alnSEXP);
    Rcpp::traits::input_parameter< std::string& >::type seq(seqSEXP);
    rcpp_result_gen = Rcpp::wrap(map_alignment_sequence(aln, seq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_synal_calc_identity", (DL_FUNC) &_synal_calc_identity, 2},
    {"_synal_map_alignment_sequence", (DL_FUNC) &_synal_map_alignment_sequence, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_synal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
