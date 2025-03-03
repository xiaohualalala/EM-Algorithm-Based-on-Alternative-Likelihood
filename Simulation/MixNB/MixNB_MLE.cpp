// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

NumericVector reorder(NumericVector vec, NumericVector ord){ 
  int n(vec.size());
  NumericVector reord(n);
  for(int i=0;i<n;i++) reord[i]=vec[ord[i]-1];
  return reord;
}
arma::vec vec2arma(Rcpp::NumericVector x)
{
    return arma::vec(x.begin(), x.size(), false); // this code originate from P159 Seamless Rcpp
}
template <typename T>
Rcpp::NumericVector arma2vec(const T &x)
{
    return Rcpp::NumericVector(x.begin(), x.end());
}
mat EmisPr(NumericVector X, NumericVector eSize, NumericVector eProb)
{
    int nbT = X.size(), nbS = eProb.size();
    mat res(nbS, nbT);

    for (unsigned i = 0; i < nbS; i++)
    {
        for (unsigned j = 0; j < nbT; j++)
        {
            res(i, j) = R::dnbinom(X(j), eSize(i), eProb(i), false);
        }
    }
    return res;
}
mat digamma(mat X)
{
    mat res(X.n_rows, X.n_cols);
    for (unsigned int i = 0; i < X.n_rows; i++)
    {
        for (unsigned int j = 0; j < X.n_cols; j++)
        {
            double tmp = R::digamma(X(i, j));
            res(i, j) = tmp;
        }
    }
    return res;
}
/*-------------------------------------------------------*
 * Function: mixNB_MLE
 * Purpose: Perform maximum likelihood estimation for the mixed negative binomial distribution using an EM algorithm.
 * Input:
 *   X (NumericVector) - Sample data.
 *   nbS (int) - Number of subpopulations or components in the mixture.
 *   avg (NumericVector, optional) - Initial mean values for each component. If not provided, k-means clustering is used for initialization.
 *   var (NumericVector, optional) - Initial variance values for each component. Used alongside 'avg' if provided.
 *   iterMax (int) - Maximum number of iterations for the EM algorithm.
 *   EPSILON (double) - Convergence threshold for the EM algorithm.
 * Output:
 *   List - A list containing estimated parameters including weights, sizes, probabilities, means, and variances for each component.
 * Description:
 *   This function estimates the parameters of a mixed negative binomial distribution, which is a mixture of several negative binomial distributions.
 *   The EM algorithm is used to iteratively estimate the parameters based on the provided data.
 *   The function handles initialization and provides detailed output of the estimated parameters.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List mixNB_MLE(NumericVector X, int nbS, NumericVector avg = NumericVector(), NumericVector var = NumericVector(), int iterMax = 15000, double EPSILON = 1e-7, bool verbose = false)
{
    // Start time for performance measurement
    Rcpp::Function SysTime("Sys.time");
    NumericVector start_time = SysTime();

    int nbT = X.size(); // Number of samples
    // Initialization of parameters
    NumericVector eWeight(nbS, 1.0 / nbS);
    NumericVector esAvg(nbS);
    NumericVector esVar(nbS);

    // Initialize means and variances if not provided
    if (avg.size() == 0 || var.size() == 0)
    {
        Rcpp::Function kmeans("kmeans");
        List res_kmeans = kmeans(X, nbS); // Perform k-means clustering for initialization

        NumericVector centers = res_kmeans[1];
        NumericVector withinss = res_kmeans[3];
        NumericVector size = res_kmeans[6];

        esAvg = clone(centers);
        esVar = withinss / (size - 1);
    }
    else
    {
        esAvg = clone(avg);
        esVar = clone(var);
    }

    // Check for invalid initial values
    if (is_true(any(esAvg >= esVar)))
    {
        stop("Invalid data: Variance should be greater than mean for each cluster.");
    }

    // Calculate initial size and probability for each component
    NumericVector eSize = pow(esAvg, 2) / (esVar - esAvg);
    NumericVector eProb = esAvg / esVar;

    // Prepare matrices for the EM algorithm
    mat ones_nbT = ones(nbT, 1);
    mat ones_nbS = ones(nbS, 1);
    mat ones_nbSS = ones(nbS, nbS);

    // Iteration starts
    int iter;
    for (iter = 0; iter < iterMax; iter++)
    {
        mat prevParam = join_cols(vec2arma(esAvg), vec2arma(eSize), vec2arma(eWeight));

        // Expectation calculation
        mat Xmat = ones_nbS * vec2arma(X).t();
        mat Nmat = vec2arma(eSize) * ones_nbT.t();

        mat Hmat = EmisPr(X, eSize, eProb);                   // Emission probabilities
        mat Pmat = Hmat % (vec2arma(eWeight) * ones_nbT.t()); // Weighted probabilities
        mat Tmat = Pmat / (ones_nbSS * Pmat);                 // Normalized probabilities

        mat Dmat = (digamma(Xmat + Nmat) - digamma(Nmat)) % Nmat; // Digamma differences for each sample and component

        NumericVector beta = 1 - 1 / (1 - eProb) - 1 / log(eProb); // Beta values for each component
        mat Bmat = vec2arma(beta) * ones_nbT.t();

        // M-Step: Maximization to update parameters
        mat eWeight_mat = (Tmat / (1.0 * nbT)) * ones_nbT;

        mat denominator = (Tmat % Dmat % Bmat) * ones_nbT;
        mat numerator = (Tmat % (Xmat + Dmat % Bmat - Dmat)) * ones_nbT;
        mat eProb_mat = denominator / numerator;

        mat temp_mat = ((Tmat % Dmat) * ones_nbT) / (Tmat * ones_nbT);
        mat eSize_mat = -temp_mat / log(eProb_mat);

        mat esAvg_mat = eSize_mat % ((1 - eProb_mat) / eProb_mat);
        mat esVar_mat = esAvg_mat / eProb_mat;

        eWeight = arma2vec(eWeight_mat); // new value of weight
        eProb = arma2vec(eProb_mat);     // new value of prob
        eSize = arma2vec(eSize_mat);     // new value of size
        esAvg = arma2vec(esAvg_mat);     // new value of mean
        esVar = arma2vec(esVar_mat);     // new value of variance

        if (any(Rcpp::is_na(esAvg)) | any(Rcpp::is_na(eSize)) | any(Rcpp::is_na(eWeight)))
        {
            warning("mixNB_MLE - NA is appeared!");
            esAvg = prevParam.col(0);
            eSize = prevParam.col(1);
            eWeight = prevParam.col(2);
            break;
        }

        Rcpp::Function order("order");
        NumericVector ord = order(eWeight); // order of mean
        esAvg = reorder(esAvg, ord);      
        eSize = reorder(eSize, ord);      
        eWeight = reorder(eWeight, ord);  
        eProb = reorder(eProb, ord);      
        esVar = reorder(esVar, ord);      

        // Check for convergence
        mat newParam = join_cols(esAvg_mat, eSize_mat, eWeight_mat);
        double D = arma::accu(arma::square(newParam - prevParam));
        if (D < EPSILON)
        {
            break;
        }
    }

    // return with a list
    List result;
    List parameter;
    parameter["eWeight"] = eWeight;
    parameter["eSize"] = eSize;
    parameter["eProb"] = eProb;
    parameter["esAvg"] = esAvg;
    parameter["esVar"] = esVar;
    result["Parameter_Estimates"] = parameter;

    // Record stop time and print results if verbose
    NumericVector stop_time = SysTime();
    Rcpp::Function difftime("difftime");
    NumericVector minutes = difftime(stop_time, start_time, "mins");
    if (verbose)
    {
        Rcpp::Rcout << "Computational iterations: " << iter + 1 << "\n";
        Rcpp::Rcout << "Computational time: " << minutes << " minutes\n";
    }

    List stats;
    stats["Iterations"] = iter + 1;
    stats["ComputationTime"] = std::to_string(minutes[0]) + " minutes";
    result["Statistics"] = stats;

    // Return with a list
    return result;
}
