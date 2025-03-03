// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#include "../Samples/Generate_Samples.h"
#include "../../MLE/C/MixNB_MLE.h"

// [[Rcpp::export]]
NumericVector mixNB_density(NumericVector x, NumericVector sizes, NumericVector probs, NumericVector weights)
{
    int n = x.size();
    NumericVector density(n, 0.0); // Initialized with 0.0

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < weights.size(); ++j)
        {
            density[i] += weights[j] * R::dnbinom(x[i], sizes[j], probs[j], false);
        }
    }

    return density;
}

// [[Rcpp::export]]
List SimulateMixNB(NumericVector probs, NumericVector sizes, NumericVector weights, int sample_size)
{
    // Generate samples
    List MixNB = SamplesMixNB(probs, sizes, weights, 1, sample_size);
    List observations_list = MixNB["samples"];
    List true_labels_list = MixNB["labels"];
    NumericVector observations = observations_list[0];
    NumericVector true_labels = true_labels_list[0];

    int K = weights.size();
    int N = observations.size();
    NumericVector predicted_labels(N);
    NumericVector prediction(N);

    // Estimate parameters using MLE
    List result = mixNB_MLE(observations, K);
    List estimation = result["Parameter_Estimates"];
    NumericVector eSize = estimation["eSize"];
    NumericVector eProb = estimation["eProb"];
    NumericVector eWeight = estimation["eWeight"];
    NumericMatrix postProbs(N, K);

    List stats = result["Statistics"];
    int iterations = stats["Iterations"];

    // Predict labels
    for (int i = 0; i < N; ++i)
    {
        // Calculate the probability for each component
        NumericVector probability(K);
        NumericVector temp(K);
        for (int k = 0; k < K; ++k)
        {
            temp[k] = R::dnbinom(observations[i], eSize[k], eProb[k], false) * eWeight[k];
        }

        double total_prob = sum(temp);
        for (int k = 0; k < K; ++k)
        {
            probability[k] = temp[k] / total_prob;
            postProbs(i, k) = probability[k];
        }
        
        // Find the index of the maximum probability
        predicted_labels[i] = which_max(probability) + 1; // C++ index starts at 0, R starts at 1
    }
    prediction = postProbs(_, 0);

    // Calculate accuracy
    double accuracy = sum(true_labels == predicted_labels) / (1.0*N);

    return List::create(Named("observations") = observations,
                        Named("estimation") = estimation,
                        Named("postProbs") = postProbs,
                        Named("true_labels") = true_labels,
                        Named("predicted_labels") = predicted_labels, 
                        Named("iterations") = iterations, 
                        Named("accuracy") = accuracy);
}


// [[Rcpp::export]]
List SimulateMixNBs(NumericVector probs, NumericVector sizes, NumericVector weights, int num_samples, int sample_size){
    NumericVector accuracy(num_samples);
    NumericVector iterations(num_samples);
    List eSize(num_samples);
    List eProb(num_samples);
    List eWeight(num_samples);

    for (int i = 0; i < num_samples; ++i){
        List res = SimulateMixNB(probs, sizes, weights, sample_size);
        accuracy[i] = res["accuracy"];
        iterations[i] = res["iterations"];
        List estimation = res["estimation"];
        eSize[i] = estimation["eSize"];
        eProb[i] = estimation["eProb"];
        eWeight[i] = estimation["eWeight"];
    }
    return List::create(Named("accuracy") = accuracy, Named("iterations") = iterations, Named("eSize") = eSize, Named("eProb") = eProb, Named("eWeight") = eWeight);
}