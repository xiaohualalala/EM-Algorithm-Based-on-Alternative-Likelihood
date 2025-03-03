# 不同方法的比较
source("../MLE/R/NB_MLE.r")
source("../Simulation/Newton-Raphson.r")
source("../Simulation/Scoring.r")
source("../Simulation/EM.r")

compare_methods <- function(prob_values, size_values, starting_values = NULL, num_samples = 200, sample_size = 2000) {
    df <- data.frame()

    for (prob in prob_values) {
        for (size in size_values) {
            samples <- SamplesNB(prob, size, num_samples, sample_size)

            new_size <- numeric(num_samples)
            new_prob <- numeric(num_samples)
            new_iterations <- numeric(num_samples)
            new_duration <- numeric(num_samples)
            new_converged <- numeric(num_samples)

            scoring_size <- numeric(num_samples)
            scoring_prob <- numeric(num_samples)
            scoring_iterations <- numeric(num_samples)
            scoring_duration <- numeric(num_samples)
            scoring_converged <- numeric(num_samples)

            em_size <- numeric(num_samples) 
            em_prob <- numeric(num_samples)
            em_iterations <- numeric(num_samples)
            em_duration <- numeric(num_samples)
            em_converged <- numeric(num_samples)

            alem_size <- numeric(num_samples)
            alem_prob <- numeric(num_samples)
            alem_iterations <- numeric(num_samples)
            alem_duration <- numeric(num_samples)
            alem_converged <- numeric(num_samples)

            for (i in 1:num_samples) {
                data <- samples[[i]]
                starting_values <- c(runif(1, 0.1, 50), runif(1, 0.01, 0.99))
                new_res <- newton_raphson(data, starting_values, max_iter = 5000)
                scoring_res <- fisher_scoring(data, starting_values, max_iter = 5000)
                em_res <- EM(data, starting_values, max_iter = 5000)
                alem_res <- NB_MLE(data, iterMax = 5000, init_values = starting_values, verbose = FALSE)

                new_size[i] <- new_res$params[1]
                new_prob[i] <- new_res$params[2]
                new_iterations[i] <- new_res$iterations
                new_duration[i] <- new_res$duration
                new_converged[i] <- new_res$converged

                scoring_size[i] <- scoring_res$params[1]
                scoring_prob[i] <- scoring_res$params[2]
                scoring_iterations[i] <- scoring_res$iterations
                scoring_duration[i] <- scoring_res$duration
                scoring_converged[i] <- scoring_res$converged

                em_size[i] <- em_res$params[1]
                em_prob[i] <- em_res$params[2]
                em_iterations[i] <- em_res$iterations
                em_duration[i] <- em_res$duration
                em_converged[i] <- em_res$converged

                alem_size[i] <- unlist(alem_res$Parameter["eSize"])
                alem_prob[i] <- unlist(alem_res$Parameter["eProb"])
                alem_iterations[i] <- alem_res$iter
                alem_duration[i] <- alem_res$duration
                alem_converged[i] <- alem_res$converged
            }
            valid_new <- which(new_converged==1)
            valid_scoring <- which(scoring_converged==1)
            valid_em <- which(em_converged==1)
            valid_alem <- which(alem_converged==1)

            temp_new <- data.frame(
                "Prob" = prob,
                "Size" = size,
                "Method" = "Newton",
                "Probability" = sum(new_converged)/length(new_converged),
                "Prob Error mean" = mean(new_prob[valid_new]),
                "Prob Error dev" = sd(new_prob[valid_new]),
                "Size Error mean" = mean(new_size[valid_new]),
                "Size Error mean" = sd(new_size[valid_new]),
                "CPU time mean" = mean(new_duration[valid_new]),
                "CPU time dev" = sd(new_duration[valid_new]),
                "Iterations mean" = mean(new_iterations[valid_new]),
                "Iterations dev" = sd(new_iterations[valid_new])
            )

            temp_scoring <- data.frame(
                "Prob" = prob,
                "Size" = size,
                "Method" = "Scoring",
                "Probability" = sum(scoring_converged)/length(scoring_converged),
                "Prob Error mean" = mean(scoring_prob[valid_scoring]),
                "Prob Error dev" = sd(scoring_prob[valid_scoring]),
                "Size Error mean" = mean(scoring_size[valid_scoring]),
                "Size Error mean" = sd(scoring_size[valid_scoring]),
                "CPU time mean" = mean(scoring_duration[valid_scoring]),
                "CPU time dev" = sd(scoring_duration[valid_scoring]),
                "Iterations mean" = mean(scoring_iterations[valid_scoring]),
                "Iterations dev" = sd(scoring_iterations[valid_scoring])
            )


            temp_em <- data.frame(
                "Prob" = prob,
                "Size" = size,
                "Method" = "EM",
                "Probability" = sum(em_converged)/length(em_converged),
                "Prob Error mean" = mean(em_prob[valid_em]),
                "Prob Error dev" = sd(em_prob[valid_em]),
                "Size Error mean" = mean(em_size[valid_em]),
                "Size Error mean" = sd(em_size[valid_em]),
                "CPU time mean" = mean(em_duration[valid_em]),
                "CPU time dev" = sd(em_duration[valid_em]),                
                "Iterations mean" = mean(em_iterations[valid_em]),
                "Iterations dev" = sd(em_iterations[valid_em])

            )

            temp_alem <- data.frame(
                "Prob" = prob,
                "Size" = size,
                "Method" = "ALEM",
                "Probability" = sum(alem_converged)/length(alem_converged),
                "Prob Error mean" = mean(unlist(alem_prob[valid_alem])),
                "Prob Error dev" = sd(unlist(alem_prob[valid_alem])),
                "Size Error mean" = mean(unlist(alem_size[valid_alem])),
                "Size Error mean" = sd(unlist(alem_size[valid_alem])),
                "CPU time mean" = mean(unlist(alem_duration[valid_alem])),
                "CPU time dev" = sd(unlist(alem_duration[valid_alem])),
                "Iterations mean" = mean(unlist(alem_iterations[valid_alem])),
                "Iterations dev" = sd(unlist(alem_iterations[valid_alem]))
            )

            temp1 <- rbind(temp_new, temp_scoring)
            temp2 <- rbind(temp1, temp_em)
            temp <- rbind(temp2,temp_alem)

            df <- rbind(df, temp)
        }
    }
    return(df)
}

set.seed(234)

prob_values <- c(0.8)
size_values <- c(8)
num_samples <- 200
sample_size <- 2000
res <- compare_methods(prob_values, size_values)
table1 <- xtable(res)

prob_values <- c(0.2, 0.5, 0.8)
size_values <- c(2, 5, 8)

# test 
newton_raphson(data)
fisher_scoring(data)
EM(data)
NB_MLE(data, verbose = FALSE)
