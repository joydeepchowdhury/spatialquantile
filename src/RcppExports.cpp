// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _spatialquantile_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _spatialquantile_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _spatialquantile_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _spatialquantile_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// g_function_weighted
double g_function_weighted(arma::mat X_local, arma::vec Q_local, arma::vec weights_local, arma::vec u_local);
RcppExport SEXP _spatialquantile_g_function_weighted(SEXP X_localSEXP, SEXP Q_localSEXP, SEXP weights_localSEXP, SEXP u_localSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_local(X_localSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Q_local(Q_localSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights_local(weights_localSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u_local(u_localSEXP);
    rcpp_result_gen = Rcpp::wrap(g_function_weighted(X_local, Q_local, weights_local, u_local));
    return rcpp_result_gen;
END_RCPP
}
// spatialquantile
arma::vec spatialquantile(arma::mat Data, arma::vec Weights, int u_index, double c, arma::vec t_vector);
RcppExport SEXP _spatialquantile_spatialquantile(SEXP DataSEXP, SEXP WeightsSEXP, SEXP u_indexSEXP, SEXP cSEXP, SEXP t_vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Weights(WeightsSEXP);
    Rcpp::traits::input_parameter< int >::type u_index(u_indexSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_vector(t_vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(spatialquantile(Data, Weights, u_index, c, t_vector));
    return rcpp_result_gen;
END_RCPP
}
// wsdrank
arma::vec wsdrank(arma::mat X_to_rank, arma::mat X_data, arma::mat X_data_weights, arma::vec t_vector);
RcppExport SEXP _spatialquantile_wsdrank(SEXP X_to_rankSEXP, SEXP X_dataSEXP, SEXP X_data_weightsSEXP, SEXP t_vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_to_rank(X_to_rankSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_data(X_dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_data_weights(X_data_weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_vector(t_vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(wsdrank(X_to_rank, X_data, X_data_weights, t_vector));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatialquantile_rcpparma_hello_world", (DL_FUNC) &_spatialquantile_rcpparma_hello_world, 0},
    {"_spatialquantile_rcpparma_outerproduct", (DL_FUNC) &_spatialquantile_rcpparma_outerproduct, 1},
    {"_spatialquantile_rcpparma_innerproduct", (DL_FUNC) &_spatialquantile_rcpparma_innerproduct, 1},
    {"_spatialquantile_rcpparma_bothproducts", (DL_FUNC) &_spatialquantile_rcpparma_bothproducts, 1},
    {"_spatialquantile_g_function_weighted", (DL_FUNC) &_spatialquantile_g_function_weighted, 4},
    {"_spatialquantile_spatialquantile", (DL_FUNC) &_spatialquantile_spatialquantile, 5},
    {"_spatialquantile_wsdrank", (DL_FUNC) &_spatialquantile_wsdrank, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatialquantile(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
