#' DEsingle: Detecting differentially expressed genes from scRNA-seq data
#'
#' This function is used to detect differentially expressed genes between two specified groups of cells in a raw read counts matrix of single-cell RNA-seq (scRNA-seq) data. It takes a non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#'
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object which contains the read counts matrix. The rows of the matrix are genes and columns are samples/cells.
#' @param group A vector of factor which specifies the two groups to be compared, corresponding to the columns in the counts matrix.
#' @param parallel If FALSE (default), no parallel computation is used; if TRUE, parallel computation using \code{BiocParallel}, with argument \code{BPPARAM}.
#' @param BPPARAM An optional parameter object passed internally to \code{\link{bplapply}} when \code{parallel=TRUE}. If not specified, \code{\link{bpparam}()} (default) will be used.
#' @return
#' A data frame containing the differential expression (DE) analysis results, rows are genes and columns contain the following items:
#' \itemize{
#'   \item theta_1, theta_2, mu_1, mu_2, size_1, size_2, prob_1, prob_2: MLE of the zero-inflated negative binomial distribution's parameters of group 1 and group 2.
#'   \item total_mean_1, total_mean_2: Mean of read counts of group 1 and group 2.
#'   \item foldChange: total_mean_1/total_mean_2.
#'   \item norm_total_mean_1, norm_total_mean_2: Mean of normalized read counts of group 1 and group 2.
#'   \item norm_foldChange: norm_total_mean_1/norm_total_mean_2.
#'   \item chi2LR1: Chi-square statistic for hypothesis testing of H0.
#'   \item pvalue_LR2: P value of hypothesis testing of H20 (Used to determine the type of a DE gene).
#'   \item pvalue_LR3: P value of hypothesis testing of H30 (Used to determine the type of a DE gene).
#'   \item FDR_LR2: Adjusted P value of pvalue_LR2 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item FDR_LR3: Adjusted P value of pvalue_LR3 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item pvalue: P value of hypothesis testing of H0 (Used to determine whether a gene is a DE gene).
#'   \item pvalue.adj.FDR: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method (Used to determine whether a gene is a DE gene).
#'   \item Remark: Record of abnormal program information.
#' }
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{DEtype}}, for the classification of differentially expressed genes found by \code{\link{DEsingle}}.
#'
#' \code{\link{TestData}}, a test dataset for DEsingle.
#'
#' @examples
#' # Load test data for DEsingle
#' data(TestData)
#'
#' # Specifying the two groups to be compared
#' # The sample number in group 1 and group 2 is 50 and 100 respectively
#' group <- factor(c(rep(1, 50), rep(2, 100)))
#'
#' # Detecting the differentially expressed genes
#' results <- DEsingle(counts = counts, group = group)
#'
#' # Dividing the differentially expressed genes into 3 categories
#' results.classified <- DEtype(results = results, threshold = 0.05)
#'
#' @import stats
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom maxLik maxLik
#' @importFrom pscl zeroinfl
#' @importMethodsFrom Matrix colSums
#' @export

DEalem <- function(counts, group, parallel = FALSE, BPPARAM = bpparam()) {
    # Handle SingleCellExperiment
    if (class(counts)[1] == "SingleCellExperiment") {
        if (!require(SingleCellExperiment)) {
            stop("To use SingleCellExperiment as input, you should install the package firstly")
        }
        counts <- counts(counts)
    }

    # Invalid input control
    if (!is.matrix(counts) & !is.data.frame(counts) & class(counts)[1] != "dgCMatrix") {
        stop("Wrong data type of 'counts'")
    }
    if (sum(is.na(counts)) > 0) {
        stop("NA detected in 'counts'")
    }
    gc()
    if (sum(counts < 0) > 0) {
        stop("Negative value detected in 'counts'")
    }
    gc()
    if (all(counts == 0)) {
        stop("All elements of 'counts' are zero")
    }
    gc()
    if (any(colSums(counts) == 0)) {
        warning("Library size of zero detected in 'counts'")
    }
    gc()

    if (!is.factor(group)) {
        stop("Data type of 'group' is not factor")
    }
    if (length(levels(group)) != 2) {
        stop("Levels number of 'group' is not two")
    }
    if (table(group)[1] < 2 | table(group)[2] < 2) {
        stop("Too few samples (< 2) in a group")
    }
    if (ncol(counts) != length(group)) {
        stop("Length of 'group' must equal to column number of 'counts'")
    }

    if (!is.logical(parallel)) {
        stop("Data type of 'parallel' is not logical")
    }
    if (length(parallel) != 1) {
        stop("Length of 'parallel' is not one")
    }

    # Preprocessing
    counts <- round(as.matrix(counts))
    storage.mode(counts) <- "integer"
    if (any(rowSums(counts) == 0)) {
        message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
    }
    counts <- counts[rowSums(counts) != 0, ]
    geneNum <- nrow(counts)
    sampleNum <- ncol(counts)
    gc()

    # Normalization
    message("Normalizing the data")
    GEOmean <- rep(NA, geneNum)
    for (i in 1:geneNum)
    {
        gene_NZ <- counts[i, counts[i, ] > 0]
        GEOmean[i] <- exp(sum(log(gene_NZ), na.rm = TRUE) / length(gene_NZ))
    }
    S <- rep(NA, sampleNum)
    counts_norm <- counts
    for (j in 1:sampleNum)
    {
        sample_j <- counts[, j] / GEOmean
        S[j] <- median(sample_j[which(sample_j != 0)])
        counts_norm[, j] <- counts[, j] / S[j]
    }
    counts_norm <- ceiling(counts_norm)
    remove(GEOmean, gene_NZ, S, sample_j, i, j)
    gc()

    # Cache totalMean and foldChange for each gene
    totalMean_1 <- rowMeans(counts[row.names(counts_norm), group == levels(group)[1]])
    totalMean_2 <- rowMeans(counts[row.names(counts_norm), group == levels(group)[2]])
    foldChange <- totalMean_1 / totalMean_2
    All_Mean_FC <- cbind(totalMean_1, totalMean_2, foldChange)

    # Memory management
    remove(counts, totalMean_1, totalMean_2, foldChange)
    counts_norm <- Matrix(counts_norm, sparse = TRUE)
    gc()


    # Function of testing homogeneity of two ZINB populations
    CallDE <- function(i) {
        # Memory management
        if (i %% 100 == 0) {
            gc()
        }

        # Function input and output
        counts_1 <- counts_norm[i, group == levels(group)[1]]
        counts_2 <- counts_norm[i, group == levels(group)[2]]
        results_gene <- data.frame(row.names = row.names(counts_norm)[i], theta_1 = NA, theta_2 = NA, mu_1 = NA, mu_2 = NA, size_1 = NA, size_2 = NA, prob_1 = NA, prob_2 = NA, total_mean_1 = NA, total_mean_2 = NA, foldChange = NA, norm_total_mean_1 = NA, norm_total_mean_2 = NA, norm_foldChange = NA, chi2LR1 = NA, pvalue_LR2 = NA, pvalue_LR3 = NA, FDR_LR2 = NA, FDR_LR3 = NA, pvalue = NA, pvalue.adj.FDR = NA, Remark = NA)

        # Log likelihood functions
        logL <- function(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) {
            logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1, pstr0 = theta_1, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2, pstr0 = theta_2, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL2 <- function(param) {
            theta_resL2 <- param[1]
            size_1_resL2 <- param[2]
            prob_1_resL2 <- param[3]
            size_2_resL2 <- param[4]
            prob_2_resL2 <- param[5]
            logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL2NZ <- function(param) {
            theta_resL2 <- 0
            size_1_resL2 <- param[1]
            prob_1_resL2 <- param[2]
            size_2_resL2 <- param[3]
            prob_2_resL2 <- param[4]
            logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3 <- function(param) {
            theta_1_resL3 <- param[1]
            size_resL3 <- param[2]
            prob_resL3 <- param[3]
            theta_2_resL3 <- param[4]
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3NZ1 <- function(param) {
            theta_1_resL3 <- 0
            size_resL3 <- param[1]
            prob_resL3 <- param[2]
            theta_2_resL3 <- param[3]
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3NZ2 <- function(param) {
            theta_1_resL3 <- param[1]
            size_resL3 <- param[2]
            prob_resL3 <- param[3]
            theta_2_resL3 <- 0
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3AZ1 <- function(param) {
            theta_1_resL3 <- 1
            size_resL3 <- param[1]
            prob_resL3 <- param[2]
            theta_2_resL3 <- param[3]
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3AZ2 <- function(param) {
            theta_1_resL3 <- param[1]
            size_resL3 <- param[2]
            prob_resL3 <- param[3]
            theta_2_resL3 <- 1
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3NZ1AZ2 <- function(param) {
            theta_1_resL3 <- 0
            size_resL3 <- param[1]
            prob_resL3 <- param[2]
            theta_2_resL3 <- 1
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        logL3NZ2AZ1 <- function(param) {
            theta_1_resL3 <- 1
            size_resL3 <- param[1]
            prob_resL3 <- param[2]
            theta_2_resL3 <- 0
            logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        judgeParam <- function(param) {
            if ((param >= 0) & (param <= 1)) {
                res <- TRUE
            } else {
                res <- FALSE
            }
            res
        }

        # MLE of parameters of ZINB counts_1
        if (sum(counts_1 == 0) > 0) {
            if (sum(counts_1 == 0) == length(counts_1)) {
                theta_1 <- 1
                mu_1 <- 0
                size_1 <- 1
                prob_1 <- size_1 / (size_1 + mu_1)
            } else {
                options(show.error.messages = FALSE)
                options(warn = -1)
                zinb_try <- try(ZINB_MLE(counts_1), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(zinb_try)) {
                    options(show.error.messages = FALSE)
                    options(warn = -1)
                    zinb_try_twice <- try(ZINB_Reg(formula_count = counts_1 ~ 1, formula_zero = counts_1 ~ 1, data = as.data.frame(counts_1), verbose = FALSE), silent = TRUE)
                    options(warn = 0)
                    options(show.error.messages = TRUE)
                    if ("try-error" %in% class(zinb_try_twice)) {
                        print("MLE of ZINB failed!")
                        results_gene[1, "Remark"] <- "ZINB failed!"
                        return(results_gene)
                    } else {
                        zinb_1 <- zinb_try_twice
                        theta_1 <- plogis(zinb_1$coefficients_zero);names(theta_1) <- NULL
                        mu_1 <- exp(zinb_1$coefficients_count);names(mu_1) <- NULL
                        prob_1 <- zinb_1$prob;names(prob_1) <- NULL
                        size_1 <- prob_1 / (1 - prob_1) * mu_1;names(size_1) <- NULL
                    }
                } else {
                    zinb_1 <- zinb_try
                    theta_1 <- zinb_1$Parameter_Estimates$eZero;names(theta_1) <- NULL
                    size_1 <- zinb_1$Parameter_Estimates$eSize;names(size_1) <- NULL
                    prob_1 <- zinb_1$Parameter_Estimates$eProb;names(prob_1) <- NULL
                    mu_1 <- (1-prob_1)/prob_1*size_1;names(mu_1) <- NULL
                }
            }
        } else {
            options(show.error.messages = FALSE)
            options(warn = -1)
            nb_try <- try(NB_MLE(counts_1), silent = TRUE)
            options(warn = 0)
            options(show.error.messages = TRUE)
            if ("try-error" %in% class(nb_try)) {
                options(show.error.messages = FALSE)
                options(warn = -1)
                nb_try_twice <- try(NB_Reg(formula = counts_1 ~ 1, data = as.data.frame(counts_1), verbose = FALSE), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(nb_try_twice)) {
                    print("MLE of NB failed!")
                    results_gene[1, "Remark"] <- "NB failed!"
                    return(results_gene)
                } else {
                    nb_1 <- nb_try_twice
                    theta_1 <- 0
                    mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
                    prob_1 <- nb_1$prob;names(prob_1) <- NULL
                    size_1 <- prob_1/(1-prob_1)*mu_1;names(size_1) <- NULL
                }
            }else{
                nb_1 <- nb_try
                theta_1 <- 0
                size_1 <- nb_1$Parameter_Estimates$eSize;names(size_1) <- NULL
                prob_1 <- nb_1$Parameter_Estimates$eProb;names(prob_1) <- NULL
                mu_1 <- (1-prob_1)/prob_1*size_1;names(mu_1) <- NULL
            }
        }

        # MLE of parameters of ZINB counts_2
        if (sum(counts_2 == 0) > 0) {
            if (sum(counts_2 == 0) == length(counts_2)) {
                theta_2 <- 1
                mu_2 <- 0
                size_2 <- 1
                prob_2 <- size_2 / (size_2 + mu_2)
            } else {
                options(show.error.messages = FALSE)
                options(warn = -1)
                zinb_try <- try(ZINB_MLE(counts_2), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(zinb_try)) {
                    options(show.error.messages = FALSE)
                    options(warn = -1)
                    zinb_try_twice <- try(ZINB_Reg(formula_count = counts_2 ~ 1, formula_zero = counts_2 ~ 1, data = as.data.frame(counts_2), verbose = FALSE), silent = TRUE)
                    options(warn = 0)
                    options(show.error.messages = TRUE)
                    if ("try-error" %in% class(zinb_try_twice)) {
                        print("MLE of ZINB failed!")
                        results_gene[1, "Remark"] <- "ZINB failed!"
                        return(results_gene)
                    } else {
                        zinb_2 <- zinb_try_twice
                        theta_2 <- plogis(zinb_2$coefficients_zero);names(theta_2) <- NULL
                        mu_2 <- exp(zinb_2$coefficients_count);names(mu_2) <- NULL
                        prob_2 <- zinb_2$prob;names(prob_2) <- NULL
                        size_2 <- prob_2 / (1 - prob_2) * mu_2;names(size_2) <- NULL
                    }
                } else {
                    zinb_2 <- zinb_try
                    theta_2 <- zinb_2$Parameter_Estimates$eZero;names(theta_2) <- NULL
                    size_2 <- zinb_2$Parameter_Estimates$eSize;names(size_2) <- NULL
                    prob_2 <- zinb_2$Parameter_Estimates$eProb;names(prob_2) <- NULL
                    mu_2 <- (1-prob_2)/prob_2*size_2;names(mu_2) <- NULL
                }
            }
        } else {
            options(show.error.messages = FALSE)
            options(warn = -1)
            nb_try <- try(NB_MLE(counts_2), silent = TRUE)
            options(warn = 0)
            options(show.error.messages = TRUE)
            if ("try-error" %in% class(nb_try)) {
                options(show.error.messages = FALSE)
                options(warn = -1)
                nb_try_twice <- try(NB_Reg(formula = counts_2 ~ 1, data = as.data.frame(counts_2), verbose = FALSE), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(nb_try_twice)) {
                    print("MLE of NB failed!")
                    results_gene[1, "Remark"] <- "NB failed!"
                    return(results_gene)
                } else {
                    nb_2 <- nb_try_twice
                    theta_2 <- 0
                    mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
                    prob_2 <- nb_2$prob;names(prob_2) <- NULL
                    size_2 <- prob_2/(1-prob_2)*mu_2;names(size_2) <- NULL
                }
            }else{
                nb_2 <- nb_try
                theta_2 <- 0
                size_2 <- nb_2$Parameter_Estimates$eSize;names(size_2) <- NULL
                prob_2 <- nb_2$Parameter_Estimates$eProb;names(prob_2) <- NULL
                mu_2 <- (1-prob_2)/prob_2*size_2;names(mu_2) <- NULL
            }
        }

        # Restricted MLE under H0 (MLE of c(counts_1, counts_2))
        if (sum(c(counts_1, counts_2) == 0) > 0) {
            options(show.error.messages = FALSE)
            options(warn = -1)
            zinb_try <- try(ZINB_MLE(c(counts_1, counts_2)), silent = TRUE)
            options(warn = 0)
            options(show.error.messages = TRUE)
            if ("try-error" %in% class(zinb_try)) {
                options(show.error.messages = FALSE)
                options(warn = -1)
                zinb_try_twice <- try(ZINB_Reg(formula_count = c(counts_1, counts_2) ~ 1, formula_zero = c(counts_1, counts_2) ~ 1, data = as.data.frame(c(counts_1, counts_2)), verbose = FALSE), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(zinb_try_twice)) {
                    print("MLE of ZINB failed!")
                    results_gene[1, "Remark"] <- "ZINB failed!"
                    return(results_gene)
                } else {
                    zinb_res <- zinb_try_twice
                    theta_res <- plogis(zinb_res$coefficients_zero);names(theta_res) <- NULL
                    mu_res <- exp(zinb_res$coefficients_count);names(mu_res) <- NULL
                    prob_res <- zinb_res$prob;names(prob_res) <- NULL
                    size_res <- prob_res/(1-prob_res)*mu_res;names(size_res) <- NULL
                }
            } else {
                zinb_res <- zinb_try
                theta_res <- zinb_res$Parameter_Estimates$eZero;names(theta_res) <- NULL
                size_res <- zinb_res$Parameter_Estimates$eSize;names(size_res) <- NULL
                prob_res <- zinb_res$Parameter_Estimates$eProb;names(prob_res) <- NULL
                mu_res <- (1-prob_res)/prob_res*size_res;names(mu_res) <- NULL
            }
        } else {
            options(show.error.messages = FALSE)
            options(warn = -1)
            nb_try <- try(NB_MLE(c(counts_1, counts_2)), silent = TRUE)
            options(warn = 0)
            options(show.error.messages = TRUE)
            if ("try-error" %in% class(nb_try)) {
                options(show.error.messages = FALSE)
                options(warn = -1)
                nb_try_twice <- try(NB_Reg(formula = c(counts_1, counts_2) ~ 1, data = as.data.frame(c(counts_1, counts_2)), verbose = FALSE), silent = TRUE)
                options(warn = 0)
                options(show.error.messages = TRUE)
                if ("try-error" %in% class(nb_try_twice)) {
                    print("MLE of NB failed!")
                    results_gene[1, "Remark"] <- "NB failed!"
                    return(results_gene)
                } else {
                    nb_res <- nb_try_twice
                    theta_res <- 0
                    mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
                    prob_res <- nb_res$prob;names(prob_res) <- NULL
                    size_res <- prob_res/(1-prob_res)*mu_res;names(size_res) <- NULL
                }
            } else {
                nb_res <- nb_try
                theta_res <- 0
                size_res <- nb_res$Parameter_Estimates$eSize;names(size_res) <- NULL
                prob_res <- nb_res$Parameter_Estimates$eProb;names(prob_res) <- NULL
                mu_res <- (1-prob_res)/prob_res*size_res;names(mu_res) <- NULL
            }
        }

        # # LRT test of H0
        chi2LR1 <- 2 * (logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_res, size_res, prob_res, counts_2, theta_res, size_res, prob_res))
        pvalue <- 1 - pchisq(chi2LR1, df = 3)

        # Format output
        results_gene[1, "theta_1"] <- theta_1
        results_gene[1, "theta_2"] <- theta_2
        results_gene[1, "mu_1"] <- mu_1
        results_gene[1, "mu_2"] <- mu_2
        results_gene[1, "size_1"] <- size_1
        results_gene[1, "size_2"] <- size_2
        results_gene[1, "prob_1"] <- prob_1
        results_gene[1, "prob_2"] <- prob_2
        results_gene[1, "norm_total_mean_1"] <- mean(counts_1)
        results_gene[1, "norm_total_mean_2"] <- mean(counts_2)
        results_gene[1, "norm_foldChange"] <- results_gene[1, "norm_total_mean_1"] / results_gene[1, "norm_total_mean_2"]
        results_gene[1, "chi2LR1"] <- chi2LR1
        results_gene[1, "pvalue"] <- pvalue

        # Restricted MLE of logL2 and logL3 under H20 and H30 when pvalue <= 0.05
        if (pvalue <= 0.05) {
            if (sum(c(counts_1, counts_2) == 0) > 0) {
                options(warn = -1)
                # Restricted MLE of logL2
                A <- matrix(rbind(c(1, 0, 0, 0, 0), c(-1, 0, 0, 0, 0), c(0, 0, 1, 0, 0), c(0, 0, -1, 0, 0), c(0, 0, 0, 0, 1), c(0, 0, 0, 0, -1)), 6, 5)
                B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 0.5, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints = list(ineqA = A, ineqB = B)), silent = TRUE)
                if ("try-error" %in% class(mleL2)) {
                    mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 0, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints = list(ineqA = A, ineqB = B)), silent = TRUE)
                }
                if ("try-error" %in% class(mleL2)) {
                    mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 1, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints = list(ineqA = A, ineqB = B)), silent = TRUE)
                }
                if ("try-error" %in% class(mleL2)) {
                    A <- matrix(rbind(c(0, 1, 0, 0), c(0, -1, 0, 0), c(0, 0, 0, 1), c(0, 0, 0, -1)), 4, 4)
                    B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                    mleL2 <- maxLik(logLik = logL2NZ, start = c(size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints = list(ineqA = A, ineqB = B))
                    theta_resL2 <- 0
                    size_1_resL2 <- mleL2$estimate["size_1_resL2"];names(size_1_resL2) <- NULL
                    prob_1_resL2 <- mleL2$estimate["prob_1_resL2"];names(prob_1_resL2) <- NULL
                    size_2_resL2 <- mleL2$estimate["size_2_resL2"];names(size_2_resL2) <- NULL
                    prob_2_resL2 <- mleL2$estimate["prob_2_resL2"];names(prob_2_resL2) <- NULL
                } else {
                    theta_resL2 <- mleL2$estimate["theta_resL2"];names(theta_resL2) <- NULL
                    size_1_resL2 <- mleL2$estimate["size_1_resL2"];names(size_1_resL2) <- NULL
                    prob_1_resL2 <- mleL2$estimate["prob_1_resL2"];names(prob_1_resL2) <- NULL
                    size_2_resL2 <- mleL2$estimate["size_2_resL2"];names(size_2_resL2) <- NULL
                    prob_2_resL2 <- mleL2$estimate["prob_2_resL2"];names(prob_2_resL2) <- NULL
                }

                # Restricted MLE of logL3
                if ((sum(counts_1 == 0) > 0) & (sum(counts_2 == 0) > 0)) {
                    # logL3
                    if (sum(counts_1 == 0) == length(counts_1)) {
                        A <- matrix(rbind(c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1)), 4, 3)
                        B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3AZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- 1
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
                    } else if (sum(counts_2 == 0) == length(counts_2)) {
                        A <- matrix(rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 0, 1), c(0, 0, -1)), 4, 3)
                        B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3AZ2, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- 1
                    } else {
                        A <- matrix(rbind(c(1, 0, 0, 0), c(-1, 0, 0, 0), c(0, 0, 1, 0), c(0, 0, -1, 0), c(0, 0, 0, 1), c(0, 0, 0, -1)), 6, 4)
                        B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
                    }
                } else if (sum(counts_1 == 0) == 0) {
                    # logL3
                    if (sum(counts_2 == 0) == length(counts_2)) {
                        A <- matrix(rbind(c(0, 1), c(0, -1)), 2, 2)
                        B <- c(1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3NZ1AZ2, start = c(size_resL3 = 1, prob_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- 0
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- 1
                    } else {
                        A <- matrix(rbind(c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1)), 4, 3)
                        B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3NZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- 0
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
                    }
                } else if (sum(counts_2 == 0) == 0) {
                    # logL3
                    if (sum(counts_1 == 0) == length(counts_1)) {
                        A <- matrix(rbind(c(0, 1), c(0, -1)), 2, 2)
                        B <- c(1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3NZ2AZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- 1
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- 0
                    } else {
                        A <- matrix(rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 0, 1), c(0, 0, -1)), 4, 3)
                        B <- c(1e-10, 1 + 1e-10, 1e-10, 1 + 1e-10)
                        mleL3 <- maxLik(logLik = logL3NZ2, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5), constraints = list(ineqA = A, ineqB = B))
                        theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
                        size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
                        prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
                        theta_2_resL3 <- 0
                    }
                }
                options(warn = 0)
            } else {
                # Restricted MLE of logL2
                theta_resL2 <- 0
                size_1_resL2 <- size_1
                prob_1_resL2 <- prob_1
                size_2_resL2 <- size_2
                prob_2_resL2 <- prob_2

                # Restricted MLE of logL3
                theta_1_resL3 <- 0
                size_resL3 <- size_res
                prob_resL3 <- prob_res
                theta_2_resL3 <- 0
            }

            # Judge parameters
            if (!(judgeParam(theta_resL2) & judgeParam(prob_1_resL2) & judgeParam(prob_2_resL2))) {
                results_gene[1, "Remark"] <- "logL2 failed!"
            }
            if (!(judgeParam(theta_1_resL3) & judgeParam(theta_2_resL3) & judgeParam(prob_resL3))) {
                results_gene[1, "Remark"] <- "logL3 failed!"
            }

            # LRT test of H20 and H30
            chi2LR2 <- 2 * (logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_resL2, size_1_resL2, prob_1_resL2, counts_2, theta_resL2, size_2_resL2, prob_2_resL2))
            pvalue_LR2 <- 1 - pchisq(chi2LR2, df = 1)
            chi2LR3 <- 2 * (logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_1_resL3, size_resL3, prob_resL3, counts_2, theta_2_resL3, size_resL3, prob_resL3))
            pvalue_LR3 <- 1 - pchisq(chi2LR3, df = 2)

            # Format output
            results_gene[1, "pvalue_LR2"] <- pvalue_LR2
            results_gene[1, "pvalue_LR3"] <- pvalue_LR3
        }

        # Return results_gene
        return(results_gene)
    }


    # Call DEG gene by gene
    if (!parallel) {
        results <- matrix(data = NA, nrow = geneNum, ncol = 22, dimnames = list(row.names(counts_norm), c("theta_1", "theta_2", "mu_1", "mu_2", "size_1", "size_2", "prob_1", "prob_2", "total_mean_1", "total_mean_2", "foldChange", "norm_total_mean_1", "norm_total_mean_2", "norm_foldChange", "chi2LR1", "pvalue_LR2", "pvalue_LR3", "FDR_LR2", "FDR_LR3", "pvalue", "pvalue.adj.FDR", "Remark")))
        results <- as.data.frame(results)
        for (i in 1:geneNum) {
            cat("\r", paste0("DEalem is analyzing ", i, " of ", geneNum, " expressed genes"))
            results[i, ] <- CallDE(i)
        }
    } else {
        message("DEalem is analyzing ", geneNum, " expressed genes in parallel")
        results <- do.call(rbind, bplapply(1:geneNum, CallDE, BPPARAM = BPPARAM))
    }

    # Format output results
    results[, c("total_mean_1", "total_mean_2", "foldChange")] <- All_Mean_FC
    results[, "FDR_LR2"] <- p.adjust(results[, "pvalue_LR2"], method = "bonferroni")
    results[, "FDR_LR3"] <- p.adjust(results[, "pvalue_LR3"], method = "bonferroni")
    results[, "pvalue.adj.FDR"] <- p.adjust(results[, "pvalue"], method = "bonferroni")
    results <- results[order(results[, "chi2LR1"], decreasing = TRUE), ]

    # Abnormity control
    if (exists("lastFuncGrad") & exists("lastFuncParam")) {
        remove(lastFuncGrad, lastFuncParam, envir = .GlobalEnv)
    }
    if (sum(!is.na(results[, "Remark"])) != 0) {
        cat(paste0("\n\n ", sum(!is.na(results[, "Remark"])), " gene failed.\n\n"))
    }

    return(results)
}
