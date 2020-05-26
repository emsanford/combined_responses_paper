// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// pwm_euclidean
double pwm_euclidean(arma::mat mat1, arma::mat mat2);
RcppExport SEXP _chromVAR_pwm_euclidean(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(pwm_euclidean(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// pwm_dist_single
arma::vec pwm_dist_single(arma::mat mat1, arma::mat mat2, arma::uword min_overlap);
RcppExport SEXP _chromVAR_pwm_dist_single(SEXP mat1SEXP, SEXP mat2SEXP, SEXP min_overlapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat2(mat2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type min_overlap(min_overlapSEXP);
    rcpp_result_gen = Rcpp::wrap(pwm_dist_single(mat1, mat2, min_overlap));
    return rcpp_result_gen;
END_RCPP
}
// compute_pwm_dist
List compute_pwm_dist(List pwms, arma::uword min_overlap);
RcppExport SEXP _chromVAR_compute_pwm_dist(SEXP pwmsSEXP, SEXP min_overlapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pwms(pwmsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type min_overlap(min_overlapSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pwm_dist(pwms, min_overlap));
    return rcpp_result_gen;
END_RCPP
}
// compute_pwm_dist2
List compute_pwm_dist2(List pwms, List pwms2, arma::uword min_overlap);
RcppExport SEXP _chromVAR_compute_pwm_dist2(SEXP pwmsSEXP, SEXP pwms2SEXP, SEXP min_overlapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pwms(pwmsSEXP);
    Rcpp::traits::input_parameter< List >::type pwms2(pwms2SEXP);
    Rcpp::traits::input_parameter< arma::uword >::type min_overlap(min_overlapSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pwm_dist2(pwms, pwms2, min_overlap));
    return rcpp_result_gen;
END_RCPP
}
// row_sds
NumericVector row_sds(arma::mat& X, bool na_rm);
RcppExport SEXP _chromVAR_row_sds(SEXP XSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(row_sds(X, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// row_sds_perm
NumericVector row_sds_perm(arma::mat& X, bool na_rm);
RcppExport SEXP _chromVAR_row_sds_perm(SEXP XSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(row_sds_perm(X, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// ProbSampleReplace
arma::urowvec ProbSampleReplace(int size, arma::vec prob);
RcppExport SEXP _chromVAR_ProbSampleReplace(SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(ProbSampleReplace(size, prob));
    return rcpp_result_gen;
END_RCPP
}
// bg_sample_helper
arma::umat bg_sample_helper(arma::uvec bin_membership, arma::mat bin_p, arma::vec bin_density, arma::uword niterations);
RcppExport SEXP _chromVAR_bg_sample_helper(SEXP bin_membershipSEXP, SEXP bin_pSEXP, SEXP bin_densitySEXP, SEXP niterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type bin_membership(bin_membershipSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bin_p(bin_pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bin_density(bin_densitySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type niterations(niterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(bg_sample_helper(bin_membership, bin_p, bin_density, niterations));
    return rcpp_result_gen;
END_RCPP
}
// euc_dist
arma::mat euc_dist(arma::mat x);
RcppExport SEXP _chromVAR_euc_dist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(euc_dist(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chromVAR_pwm_euclidean", (DL_FUNC) &_chromVAR_pwm_euclidean, 2},
    {"_chromVAR_pwm_dist_single", (DL_FUNC) &_chromVAR_pwm_dist_single, 3},
    {"_chromVAR_compute_pwm_dist", (DL_FUNC) &_chromVAR_compute_pwm_dist, 2},
    {"_chromVAR_compute_pwm_dist2", (DL_FUNC) &_chromVAR_compute_pwm_dist2, 3},
    {"_chromVAR_row_sds", (DL_FUNC) &_chromVAR_row_sds, 2},
    {"_chromVAR_row_sds_perm", (DL_FUNC) &_chromVAR_row_sds_perm, 2},
    {"_chromVAR_ProbSampleReplace", (DL_FUNC) &_chromVAR_ProbSampleReplace, 2},
    {"_chromVAR_bg_sample_helper", (DL_FUNC) &_chromVAR_bg_sample_helper, 4},
    {"_chromVAR_euc_dist", (DL_FUNC) &_chromVAR_euc_dist, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_chromVAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
