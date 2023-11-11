//////////////////////////////////////////////////////////////////////////
//
// * FILE: calipeRmatching.cpp
// * DESCRIPTION:
//    Caliper matching library. 
//    Interface, to be integrated in an R package.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 17/oct/2023
//  * NOTE: To add new propensity score models, follow the instructions
//          indicated by the comment `EDIT_TO_ADD`. In particular,
//          (0) modify `propscore.c` and `propscore.h` as described
//              therein; and
//          (1) adjust `cm_number_of_modeltypes` in `cm.c`; and
//          (2) add the name of the model to the
//          `CM_FOREACH_MODEL` definition in `cm.c`; and
//          (3) add the new model to the g/ginv/gderiv/ginvderiv arrays in `cm.c`.
//          (4) add the new model in this file where indicated by EDIT_TO_ADD.
//
//////////////////////////////////////////////////////////////////////////
// #include <RcppGSL.h>
#include <Rcpp.h>
#include <thread>

#include "cm.h"


// set the default number of threads to be used; if estimate is 0, default to using 1. Can be overwritten by user using `cpp_cm_set_number_of_threads()`.
const int estimated_number_of_threads = std::thread::hardware_concurrency();
int NUM_THREADS = estimated_number_of_threads > 0 ? estimated_number_of_threads : 1;


//' Sets the number or parallel threads for calipeRmatching library.
//'
//' Only for internal use. Users should use `cm_set_number_of_threads()` that has the same funcionality.
//' Sets the global variable NUM_THREADS to `nr_threads` for parallel execution.
//' 
//' @param nr_threads the desired number of threads.
//' @return Zero if NUM_THREADS sets successfully. Otherwise, the program stops.
//' 
//' @author Mate Kormos
//' @references [KPV23] Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373
//'
// [[Rcpp::export]]
int cpp_cm_set_number_of_threads(int nr_threads){
    if (nr_threads <= 0){
        Rcpp::stop("%s: line %d: `nr_threads` must be strictly positive; now `nr_threads` = %d.\n", __FILE__, __LINE__, nr_threads);
    }
    NUM_THREADS = nr_threads;
    return 0;
}


//' Test suite for calipeRmatching library.
//'
//' Only for internal use. Users should use `test_cm()` that has the same funcionality.
//' 
//' @return Zero if all tests passed successfully. Otherwise, the program stops.
//' 
//' @author Mate Kormos
//' @references [KPV23] Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373
//'
// [[Rcpp::export]]
int cpp_test_cm(){
    // const Rcpp::Environment global_r_env = Rcpp::Environment::global_env();
    // const int estimated_number_of_threads = global_r_env["N_THREADS"];
    // Rprintf("You have %d processors.\n", std::thread::hardware_concurrency());
    // Rprintf("You have %d processors.\n", NUM_THREADS);
    test_cm();
    return 0;
}


//' Caliper matching estimator for known propensity score cpp interface.
//'
//' Only for internal use. Users should use `cm_cm_known_propscore()` that has the same functionality with error checking.
//' 
//' @param y `n`-long \code{numeric} \code{vector} of observations of the outcome variable.
//' @param d `n`-long \code{numeric} \code{vector} of observations of the treatment indicator
//'  variable. Each entry can be either 0 for control or 1 for treatment status.
//' @param propscore `n`-long \code{numeric} \code{vector} of observations of the propensity score.
//' @param delta the caliper, a positive \code{numeric} variable. If zero is passed, then the default 
//' data-driven value is used (recommended). If a positive value is passed, that is used instead.
//' @param estimate_variance \code{logical} variable. If \code{TRUE}, the variances of the estimators are estimated.
//' If \code{FALSE}, no variance estimates are provided; this is useful to speed up execution.
//' @param beta a positive \code{numeric} variable, negative-exponent of bandwidth in nonparametric variance estimation. 
//' If zero  is passed, it is the default value (recommended), otherwise it is equal to the passed value.
//'@param alpha a positive \code{numeric} variable, negative-exponent of truncation sequence in nonparametric variance estimation. 
//' If zero  is passed, it is the default value (recommended), otherwise it is equal to the passed value. 
//' Must be strictly smaller than `beta`.
//'@param kappa_a a positive \code{numeric} variable, scale of truncation sequence in nonparametric variance estimation. 
//' If zero  is passed, it is the default value, otherwise it is equal to the passed value.
//'@param kappa_gamma a positive \code{numeric} variable, scale of bandwidth sequence in nonparametric variance estimation. 
//' If zero is passed, it is the default value, otherwise it is equal to the passed value.
//' 
//'@return A \code{List}, results from caliper matching estimation with elements:
//' \itemize{
//' \item point_estimates: a \code{List} with elements \itemize{ 
//'      \item ate_hat: estimated ATE; \eqn{\hat{\tau}_{\pi}} in [KPV23].
//'      \item att_hat: estimated ATT; \eqn{\hat{\tau}_{t,\pi}} in [KPV23].
//'      \item var_hat_ate: estimated asymptotic variance of `ate_hat`; \eqn{ \hat{V}_{\tau} + \hat{V}_{\sigma,\pi} } in [KPV23].
//'      \item var_hat_att: estimated asymptotic variance of `att_hat`; \eqn{ \hat{V}_{\tau_{t}} + \hat{V}_{t,\sigma,\pi} } in [KPV23].
//'      \item var_hat_component_tau_ate: the estimate \eqn{ \hat{V}_{\tau} } in [KPV23].
//'      \item var_hat_component_tau_att: the estimate \eqn{ \hat{V}_{\tau_{t}} } in [KPV23].
//'      \item var_hat_component_sigmapi_ate: the estimate \eqn{ \hat{V}_{\sigma,\pi} } in [KPV23].
//'      \item var_hat_component_sigmapi_att: the estimate \eqn{ \hat{V}_{t,\sigma,\pi} } in [KPV23].
//'      }
//' \item delta: caliper that is actually used. If zero was passed, then it is equal to the default data-driven value. If a positive value was passed, then it is equal to that instead.
//' \item number_of_matches: `n`-long vector, the number of matches for each unit.
//' \item matched_pairs: a 2-by-(number of matched pairs) matrix. A column with value (i,j) indicates that units i and j are matched. Indices start from 1.  
//' \item var_estimation_settings: a \code{List} with elements \itemize{
//'      \item estimate_variance: whether variance estimation is performed (`estimate_variance=TRUE`), or not (`estimate_variance=FALSE`).
//'      \item a_n: value of truncation sequence in nonparametric variance estimation. Equal to `kappa_a` * `n` ^ `alpha`.
//'      \item gamma_n:  value of bandwidth in nonparametric variance estimation. Equal to `kappa_gamma` * `n` ^ `beta`.
//'      \item beta: negative-exponent of bandwidth in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item alpha: negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item kappa_gamma: scale of bandwidth sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item kappa_a: scale of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item propscore_min: the smallest propensity score value in `propscore`.
//'      \item propscore_max: the largest propensity score value in `propscore`.
//'      \item truncation_low: lower threshold for variance estimation = `propscore_min` + `a_n`.
//'      \item truncation_high: higher threshold for variance estimation = `propscore_max` - `a_n`.
//'      }
//' }
//' 
//'@section Warning: 
//' 
//' It is advised to look at the variance components in the returned results. 
//' Undesirably, negative values may occur; this may indicate too strong truncation at the bounderies.
//' To resolve this, lower `kappa_a` to mitigate truncation. 
//' Alternatively, try adjusting `kappa_gamma`.
//' 
//' @seealso \link{cm_cm_estimated_propscore}
//' @author Mate Kormos
//' @references [KPV23] Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373
//'
// [[Rcpp::export]]
Rcpp::List cpp_cm_cm_known_propscore(const Rcpp::NumericVector y, // outcome variables
                                     const Rcpp::NumericVector d, // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise   
                                     const Rcpp::NumericVector propscore, // vector of propensity score values
                                     double delta, // caliper
                                     bool estimate_variance,  // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates are not required.
                                     double beta,    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
                                     double alpha,   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
                                     double kappa_a, // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
                                     double kappa_gamma // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
    ){
// Rcpp::List cpp_cm_cm_known_propscore(const RcppGSL::vector<double> &y, // outcome variables
//                                      const RcppGSL::vector<int> &d, // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise   
//                                      const RcppGSL::vector<double> &propscore, // vector of propensity score values
//                                      double delta, // caliper
//                                      double beta,    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
//                                      double alpha,   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
//                                      double kappa_a, // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
//                                      double kappa_gamma // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
//     ){
    
    // process inputs
    // int n = y->size; // if `y` is RcppGSL::vector<double> type
    if (!((y.size() == d.size()) && (d.size() == propscore.size()))){  // sanity check
        Rprintf("%s: line %d: input check failed: `y`, `d`, `propscore` must have same length; now they have `y`.size = %td, `d`.size = %td, `propscore`.size = %td.\n", __FILE__, __LINE__, y.size(), d.size(), propscore.size());
        Rcpp::stop("");
    };
    int n = y.size();  // if `y` is Rcpp::NumericVector type
    vector *_y = vector_alloc(n);
    vector_short *_d = vector_short_alloc(n);
    vector *_propscore = vector_alloc(n);
    // --- copy data from Rcpp types to vector types to be handled by cm.c
    for (int i=0; i<n; i++){
        vector_set(_y, i, y[i]);
        vector_short_set(_d, i, d[i]);
        vector_set(_propscore, i, propscore[i]);
    }
    // --- estimate variance
    int _estimate_variance = 0;
    if (estimate_variance){
        _estimate_variance = 1;
    }

    // setup model
    CMModelKnownPropscore *cm_model_known_propscore = (CMModelKnownPropscore *) malloc(sizeof(CMModelKnownPropscore));
    cm_model_known_propscore->y = _y;
    cm_model_known_propscore->d = _d;
    cm_model_known_propscore->propscore = _propscore;
    cm_model_known_propscore->delta = delta;
    cm_model_known_propscore->estimate_variance = _estimate_variance;
    cm_model_known_propscore->beta = beta;
    cm_model_known_propscore->alpha = alpha;
    cm_model_known_propscore->kappa_gamma = kappa_gamma;
    cm_model_known_propscore->kappa_a = kappa_a;
    cm_initialise_known_propscore(cm_model_known_propscore);
    
    // estimation
    CMResults *cm_results = (CMResults *) malloc(sizeof(CMResults)); // to store estimation results
    cm_results = cm_cm_known_propscore(cm_model_known_propscore);
    
    // process results
    // --- number of matches
    Rcpp::IntegerVector number_of_matches(n);
    for (int i=0; i<n; i++){
        number_of_matches[i] = vector_int_get(cm_results->number_of_matches, i);
    }
    // --- match indices as a 2-by-(number of matched pairs) matrix; 
    // R is column-major order so flat matrix is better provided the two coordinates of a match are accessed after each other
    int number_of_matched_pairs = vector_int_sum(cm_results->number_of_matches) / 2; // every match is counted twice due to symmetry of caliper matching
    Rcpp::IntegerMatrix matched_pairs(2, number_of_matched_pairs);
    int  matched_pair_counter = 0;
    int current_matchindex_i;
    for (int i=0; i<n; i++){
        for (int j=0; j<cm_results->match_indices[i]->usage; j++){  // traverse matches of unit i
            current_matchindex_i = dynamicarray_int_get(cm_results->match_indices[i], j); // c index starting from 0
            // only add matched (i,j(i)) if it wasn't added before as (j(i),i)
            if (current_matchindex_i > i){  // add and only add (i, j(i)) if j(i) that is larger than i; if it is smaller then has been added already as (j(i), i).
                    matched_pairs(0, matched_pair_counter) = i;
                    matched_pairs(1, matched_pair_counter) = current_matchindex_i; 
                    matched_pair_counter++;
            }   
        }
    }
    matched_pairs = matched_pairs + 1; // +1 because R indexes from 1



    // cleanup
    vector_free(_y);
    vector_short_free(_d);
    vector_free(_propscore);
    
    // results as list
    return Rcpp::List::create(
        // outputs
        Rcpp::Named("point_estimates") = Rcpp::List::create(
            Rcpp::Named("ate_hat") = cm_results->ate_hat,
            Rcpp::Named("att_hat") = cm_results->att_hat,
            Rcpp::Named("var_hat_ate") = cm_results->var_hat_ate,
            Rcpp::Named("var_hat_att") = cm_results->var_hat_att,
            Rcpp::Named("var_hat_component_tau_ate") = cm_results->var_hat_component_tau_ate,  // estimated V_tau
            Rcpp::Named("var_hat_component_tau_att") = cm_results->var_hat_component_tau_att,  // estimated V_tau
            Rcpp::Named("var_hat_component_sigmapi_ate") = cm_results->var_hat_component_sigmapi_ate,  // estimated V_sigmapi
            Rcpp::Named("var_hat_component_sigmapi_att") = cm_results->var_hat_component_sigmapi_att  // estimated V_sigmapi
            ),
        Rcpp::Named("delta") = cm_results->delta, // caliper. If the passed value was strictly positive, it is equal to it, if the passed value was zero, then it is calculated from data.
        Rcpp::Named("number_of_matches") = number_of_matches, // number of matches for each unit
        Rcpp::Named("matched_pairs") = matched_pairs, // a 2-by-(number of matched pairs) matrix. A column with value (i,j) indicates that units i and j are matched. Indices start from 1.  
        // fields related to nonparametric variance estimation
        Rcpp::Named("var_estimation_settings") = Rcpp::List::create(
            Rcpp::Named("estimate_variance") = cm_results->estimate_variance,
            Rcpp::Named("a_n") = cm_results->a_n, // value of truncation sequence in nonparametric variance estimation. Equal to `kappa_a` * `n` ^ `alpha`.
            Rcpp::Named("gamma_n") = cm_results->gamma_n, // value of bandwidth in nonparametric variance estimation. Equal to `kappa_gamma` * `n` ^ `beta`.
            Rcpp::Named("beta") = cm_results->beta,  // negative-exponent of bandwidth in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("alpha") = cm_results->alpha,  // negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("kappa_gamma") = cm_results->kappa_gamma,  // scale of bandwidth sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("kappa_a") = cm_results->kappa_a,  // scale of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            // fields computed from data
            Rcpp::Named("propscore_min") = cm_results->propscore_min, // the smallest propensity score value
            Rcpp::Named("propscore_max") = cm_results->propscore_max, // the largest propensity score value
            Rcpp::Named("truncation_low") = cm_results->truncation_low, // lower threshold for variance estimation = `propscore_min`+`a_n`
            Rcpp::Named("truncation_high") = cm_results->truncation_high // higher threshold for variance estimation = `propscore_max` - `a_n`
        )
    );
}



//' Caliper matching estimator for estimated propensity score cpp interface.
//'
//' Only for internal use. Users should use `cm_cm_estimated_propscore()` that has the same functionality with error checking.
//' @param y `n`-long \code{numeric} \code{vector} of observations of the outcome variable.
//' @param d `n`-long \code{numeric} \code{vector} of observations of the treatment indicator
//'  variable. Each entry can be either 0 for control or 1 for treatment status.
//' @param x `n`-by-`k` \code{numeric} \code{array} of `n` observations of the `k` covariates.
//' @param modeltype type of the propensity score model. Currently only "logit" or "probit" are supported.
//' @param theta estimated `k+1`-long coefficient vector of the propensity score model. It is assumed to be statistically independent of (`y`,`d,`x`).
//' The first entry corresponds to the intercept.
//' @param delta the caliper, a positive \code{numeric} variable. If zero is passed, then the default 
//' data-driven value is used (recommended). If a positive value is passed, that is used instead.
//' @param estimate_variance \code{logical} variable. If \code{TRUE}, the variances of the estimators are estimated.
//' If \code{FALSE}, no variance estimates are provided; this is useful to speed up execution.
//' @param beta a positive \code{numeric} variable, negative-exponent of bandwidth in nonparametric variance estimation. 
//' If zero  is passed, it is the default value (recommended), otherwise it is equal to the passed value.
//'@param alpha a positive \code{numeric} variable, negative-exponent of truncation sequence in nonparametric variance estimation. 
//' If zero  is passed, it is the default value (recommended), otherwise it is equal to the passed value.
//' Must be strictly smaller than `beta`.
//'@param kappa_a a positive \code{numeric} variable, scale of truncation sequence in nonparametric variance estimation. 
//' If zero  is passed, it is the default value, otherwise it is equal to the passed value.
//'@param kappa_gamma a positive \code{numeric} variable, scale of bandwidth sequence in nonparametric variance estimation. 
//' If zero is passed, it is the default value, otherwise it is equal to the passed value.
//'@param kappa_gamma_derivative scale parameter of bandwidth in nonparametric variance estimation estimating derivatives w.r.t. \eqn{\theta}, the propensity score parameter. 
//' If zero is passed, a default value is used.
//' 
//'@return A \code{List}, results from caliper matching estimation with elements:
//' \itemize{
//' \item point_estimates: a \code{List} with elements \itemize{ 
//'      \item ate_hat: estimated ATE; \eqn{\hat{\tau}_{\pi}} in [KPV23].
//'      \item att_hat: estimated ATT; \eqn{\hat{\tau}_{t,\pi}} in [KPV23].
//'      \item var_hat_ate: estimated asymptotic variance of `ate_hat`; \eqn{ \hat{V}_{\tau} + \hat{V}_{\sigma,\pi} } in [KPV23].
//'      \item var_hat_att: estimated asymptotic variance of `att_hat`; \eqn{ \hat{V}_{\tau_{t}} + \hat{V}_{t,\sigma,\pi} } in [KPV23].
//'      \item var_hat_component_tau_ate: the estimate \eqn{ \hat{V}_{\tau} } in [KPV23].
//'      \item var_hat_component_tau_att: the estimate \eqn{ \hat{V}_{\tau_{t}} } in [KPV23].
//'      \item var_hat_component_sigmapi_ate: the estimate \eqn{ \hat{V}_{\sigma,\pi} } in [KPV23].
//'      \item var_hat_component_sigmapi_att: the estimate \eqn{ \hat{V}_{t,\sigma,\pi} } in [KPV23].
//'      \item var_hat_component_estpi_ate: estimated variance component deriving from the estimation of propensity score; the third term in Equation (14) in [KPV23].
//'      \item var_hat_component_estpi_att: estimated variance component deriving from the estimation of propensity score; the third term in Equation (15) in [KPV23].
//' }
//' \item delta: caliper that is actually used. If zero was passed, then it is equal to the default data-driven value. If a positive value was passed, then it is equal to that instead.
//' \item number_of_matches: `n`-long vector, the number of matches for each unit.
//' \item matched_pairs: a 2-by-(number of matched pairs) matrix. A column with value (i,j) indicates that units i and j are matched. Indices start from 1.
//' \item var_estimation_settings: a \code{List} with elements \itemize{
//'      \item estimate_variance: whether variance estimation is performed (`estimate_variance=TRUE`), or not (`estimate_variance=FALSE`).
//'      \item a_n: value of truncation sequence in nonparametric variance estimation. Equal to `kappa_a` * `n` ^ `alpha`.
//'      \item gamma_n:  value of bandwidth in nonparametric variance estimation. Equal to `kappa_gamma` * `n` ^ `beta`.
//'      \item beta: negative-exponent of bandwidth in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item alpha: negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item kappa_gamma: scale of bandwidth sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item kappa_a: scale of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
//'      \item kappa_gamma_derivative: negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.     
//'      \item propscore_min: the smallest propensity score value in `propscore`.
//'      \item propscore_max: the largest propensity score value in `propscore`.
//'      \item truncation_low: lower threshold for variance estimation = `propscore_min` + `a_n`.
//'      \item truncation_high: higher threshold for variance estimation = `propscore_max` - `a_n`.
//'      }
//' }
//' 
//'@section Warning: 
//' 
//' It is advised to look at the variance components in the returned results. 
//' Undesirably, negative values may occur; this may indicate too strong truncation at the bounderies.
//' To resolve this, lower `kappa_a` to mitigate truncation. 
//' Alternatively, try adjusting `kappa_gamma` and/or `kappa_gamma_derivative`.
//' 
//' 
//' @seealso \link{cm_cm_known_propscore}
//' @author Mate Kormos
//' @references [KPV23] Kormos, V. d. Pas, V. d. Vaart (2023): Asymptotics of Caliper Matching Estimators for Average Treatment Effects, https://arxiv.org/abs/2304.08373
//' 
//'
// [[Rcpp::export]]
Rcpp::List cpp_cm_cm(const Rcpp::NumericVector y, // outcome variable.
                     const Rcpp::NumericVector d, // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise.
                     const Rcpp::NumericMatrix x, // `n`-by-`k` matrix of covariates.
                     std::string modeltype, // propensity score modeltype
                     const Rcpp::NumericVector theta, // the estimated propensity score parameter.
                     double delta, // caliper.
                     bool estimate_variance,  // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates are not required.
                     double beta,    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
                     double alpha,   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
                     double kappa_a, // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
                     double kappa_gamma, // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
                     double kappa_gamma_derivative // scale parameter of bandwidth in nonparametric variance estimation estimating derivatives w.r.t. theta. If zero is passed, a default value is used. 
    ){   
    // process inputs
    // --- data dimensions
    if (!((y.size() == d.size()) && (d.size() == x.nrow()))){  // sanity check
        Rprintf("%s: line %d: input check failed: `y`, `d` must have same length equal to the number of rows in `x`; now `y`.size = %td, `d`.size = %td, `x`.nrow = %td.\n", __FILE__, __LINE__, y.size(), d.size(), x.nrow());
        Rcpp::stop("");
    };
    if (!( (x.ncol() + 1) == theta.size())){
        Rprintf("%s: line %d: input check failed: `x` must have exactly one more column than length of `theta`; now `x`.nrow = %td, `theta`.size = %td.\n", __FILE__, __LINE__, x.nrow(), theta.size());
        Rcpp::stop("");
    }
    int n = y.size();  // if `y` is Rcpp::NumericVector type
    int k = x.ncol();
    vector *_y = vector_alloc(n);
    vector_short *_d = vector_short_alloc(n);
    vector *_propscore = vector_alloc(n);
    matrix *_x = matrix_alloc(n, k);
    // --- copy data from Rcpp types to vector types to be handled by cm.c
    for (int i=0; i<n; i++){
        vector_set(_y, i, y[i]);
        vector_short_set(_d, i, d[i]);
        for (int j=0; j<k; j++){
            matrix_set(_x, i, j, x(i, j));
        }
    }
    // --- modeltype
    int *char_length = (int *) malloc(sizeof(int));
    *char_length = modeltype.length() + 1;
    char *_modeltype = (char *) malloc(*char_length);
    strncpy(_modeltype, modeltype.c_str(), *char_length);
    free(char_length);
    // --- propensity score parameter theta
    vector *_theta = vector_alloc(k + 1);
    vector_set(_theta, k, theta[0]); // intercept being the first entry in R is moved to last entry 
    for (int j=0; j<k; j++){
        vector_set(_theta, j, theta[j + 1]);
    }
    for (int j=0; j<k+1; j++){ // sanity check
        Rprintf("_theta[%d] = %f.\n", j, vector_get(_theta, j));
    }
    // --- estimate variance
    int _estimate_variance = 0;
    if (estimate_variance){
        _estimate_variance = 1;
    }

    // setup model
    CMModel *cm_model = (CMModel *) malloc(sizeof(CMModel));
    cm_model->y = _y;
    cm_model->d = _d;
    cm_model->x = _x;
    cm_model->delta = delta;
    cm_model->modeltype = _modeltype;
    cm_model->theta = _theta;
    cm_model->estimate_variance = _estimate_variance;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    cm_initialise(cm_model);
    
    // estimation
    CMResults *cm_results = (CMResults *) malloc(sizeof(CMResults)); // to store estimation results
    cm_results = cm_cm(cm_model);
    
    // process results
    // --- number of matches
    Rcpp::IntegerVector number_of_matches(n);
    for (int i=0; i<n; i++){
        number_of_matches[i] = vector_int_get(cm_results->number_of_matches, i);
    }
    // --- match indices as a 2-by-(number of matched pairs) matrix; 
    // R is column-major order so flat matrix is better provided the two coordinates of a match are accessed after each other
    int number_of_matched_pairs = vector_int_sum(cm_results->number_of_matches) / 2; // every match is counted twice due to symmetry of caliper matching
    Rcpp::IntegerMatrix matched_pairs(2, number_of_matched_pairs);
    int matched_pair_counter = 0;
    int current_matchindex_i;
    for (int i=0; i<n; i++){
        for (int j=0; j<cm_results->match_indices[i]->usage; j++){  // traverse matches of unit i
            current_matchindex_i = dynamicarray_int_get(cm_results->match_indices[i], j); // c index starting from 0
            // only add matched (i,j(i)) if it wasn't added before as (j(i),i)
            if (current_matchindex_i > i){  // add and only add (i, j(i)) if j(i) that is larger than i; if it is smaller then has been added already as (j(i), i).
                    matched_pairs(0,  matched_pair_counter) = i;
                    matched_pairs(1,  matched_pair_counter) = current_matchindex_i; 
                    matched_pair_counter++;
            }   
        }
    }
    matched_pairs = matched_pairs + 1; // +1 because R indexes from 1



    // cleanup
    vector_free(_y);
    vector_short_free(_d);
    matrix_free(_x);
    
    // results as list
    return Rcpp::List::create(
        // outputs
        Rcpp::Named("point_estimates") = Rcpp::List::create(
            Rcpp::Named("ate_hat") = cm_results->ate_hat,
            Rcpp::Named("att_hat") = cm_results->att_hat,
            Rcpp::Named("var_hat_ate") = cm_results->var_hat_ate,
            Rcpp::Named("var_hat_att") = cm_results->var_hat_att,
            Rcpp::Named("var_hat_component_tau_ate") = cm_results->var_hat_component_tau_ate,  // estimated V_tau
            Rcpp::Named("var_hat_component_tau_att") = cm_results->var_hat_component_tau_att,  // estimated V_tau
            Rcpp::Named("var_hat_component_sigmapi_ate") = cm_results->var_hat_component_sigmapi_ate,  // estimated V_sigmapi
            Rcpp::Named("var_hat_component_sigmapi_att") = cm_results->var_hat_component_sigmapi_att,  // estimated V_sigmapi
            Rcpp::Named("var_hat_component_estpi_ate") = cm_results->var_hat_component_estpi_ate,  // estimated variance component deriving from the estimation of prepensity score
            Rcpp::Named("var_hat_component_estpi_att") = cm_results->var_hat_component_estpi_att  // estimated variance component deriving from the estimation of prepensity score
        ),
        Rcpp::Named("delta") = cm_results->delta, // caliper. If the passed value was strictly positive, it is equal to it, if the passed value was zero, then it is calculated from data.
        Rcpp::Named("number_of_matches") = number_of_matches, // number of matches for each unit
        Rcpp::Named("matched_pairs") = matched_pairs, // a 2-by-(number of matched pairs) matrix. A column with value (i,j) indicates that units i and j are matched. Indices start from 1.  
        // fields related to nonparametric variance estimation
        Rcpp::Named("var_estimation_settings") = Rcpp::List::create(
            Rcpp::Named("estimate_variance") = cm_results->estimate_variance,
            Rcpp::Named("a_n") = cm_results->a_n, // value of truncation sequence in nonparametric variance estimation. Equal to `kappa_a` * `n` ^ `alpha`.
            Rcpp::Named("gamma_n") = cm_results->gamma_n, // value of bandwidth in nonparametric variance estimation. Equal to `kappa_gamma` * `n` ^ `beta`.
            Rcpp::Named("gamma_derivative_n") = cm_results->gamma_derivative_n, // value of bandwidth in nonparametric variance estimation used in derivative estimation w.r.t. propensity score parameters; it is zero when the propensity score is known. Equal to `kappa_gamma_derivative` * `n` ^ `beta`.
            Rcpp::Named("beta") = cm_results->beta,  // negative-exponent of bandwidth in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("alpha") = cm_results->alpha,  // negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("kappa_gamma") = cm_results->kappa_gamma,  // scale of bandwidth sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("kappa_a") = cm_results->kappa_a,  // scale of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
            Rcpp::Named("kappa_gamma_derivative") = cm_results->kappa_gamma_derivative,  // negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value (which is zero for known propensity scores), otherwise it is equal to the passed value.
            // fields computed from data
            Rcpp::Named("propscore_min") = cm_results->propscore_min, // the smallest propensity score value
            Rcpp::Named("propscore_max") = cm_results->propscore_max, // the largest propensity score value
            Rcpp::Named("truncation_low") = cm_results->truncation_low, // lower threshold for variance estimation = `propscore_min`+`a_n`
            Rcpp::Named("truncation_high") = cm_results->truncation_high // higher threshold for variance estimation = `propscore_max` - `a_n`
        )
    );
}
