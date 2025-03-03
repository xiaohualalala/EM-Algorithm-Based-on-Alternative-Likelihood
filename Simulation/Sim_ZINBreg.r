source('../Reg/ZINB_Reg.r', chdir = TRUE)
library(ggplot2)
library(gridExtra)

# ZINB Reg
generate_zinb_data <- function(n, beta, gamma, theta, X_intercept = TRUE, Z_intercept = TRUE, X_dist = rnorm, Z_dist = rnorm) {
  # n: 样本数量
  # beta: 计数部分的回归系数
  # gamma: 零膨胀部分的回归系数
  # theta: 离散参数
  # X_dist: 自变量 X 的生成分布函数
  # Z_dist: 自变量 Z 的生成分布函数
  
  # 1. 生成自变量 X 和 Z
  p_count <- length(beta) - X_intercept # 计数部分变量个数
  p_zero <- length(gamma) - Z_intercept  # 零膨胀部分变量个数
  
  X <- matrix(X_dist(n * p_count), nrow = n, ncol = p_count)
  Z <- matrix(Z_dist(n * p_zero), nrow = n, ncol = p_zero)
  if (X_intercept) {
    X <- cbind(Intercept_X = 1, X)
  }
  if (Z_intercept) {
    Z <- cbind(Intercept_Z = 1, Z)
  }

  # 2. 计算均值 μ 和零膨胀概率 π
  mu <- exp(X %*% beta)  # 计数部分均值
  pi <- 1 / (1 + exp(-Z %*% gamma))  # 零膨胀概率
  
  # 3. 模拟响应变量 Y 并记录零的来源
  Y <- numeric(n)
  source <- character(n)  # 用于记录零的来源
  
  for (i in 1:n) {
    if (runif(1) < pi[i]) {
      # 零膨胀部分：Y = 0
      Y[i] <- 0
      source[i] <- "zero_part"
    } else {
      # 计数部分：从负二项分布中生成 Y
      Y[i] <- rnbinom(1, size = theta, mu = mu[i])
      source[i] <- "count_part"
    }
  }
  
  # 返回数据框
  data <- data.frame(Y = Y, X = X[, -1, drop = FALSE], Z = Z[, -1, drop = FALSE], source = source)
  colnames(data) <- c("Y", paste0("X", 1:p_count), paste0("Z", 1:p_zero), "source")
  return(data)
}
set.seed(42)
data <- generate_zinb_data(2000, beta = c(2, -0.5), gamma = c(-1, 1.5, -0.5), theta = 2)

fit <- ZINB_Reg(
  formula_count = Y ~ X1,
  formula_zero = Y ~ Z1 + Z2,
  data = data,
  verbose = FALSE
)
fit2 <- ZINB_Reg(
  formula_count = Y ~ X1,
  formula_zero = Y ~ Z1 + Z2,
  data = data,
  type = "bWald"
)
fit$summary
fit2$summary

library(pscl)
fit3 <- zeroinfl(Y ~ X1 | Z1 + Z2, data = data, dist = "negbin")
summary(fit3)

library(glmmTMB)
fit4 <- glmmTMB(
  formula = Y ~ X1,  
  ziformula = ~ Z1 + Z2,              
  family = nbinom2,                  
  data = data
)
fit4 <- summary(model)

fit5 <- gamlssML(
  Y ~ X1 | Z1 + Z2, 
  data = data, 
  family="ZINBI")

fit5 <- gamlssML(
  formula = Y ~  X1,        # 计数部分（mu）
  nu.formula = ~ 1 + Z1 + Z2,  # 零膨胀部分（nu）
  family = ZINBI,
  data = data
)
fit4

predict2 <- function(fit, formula_count, formula_zero, data) {
    Z <- model.matrix(formula_zero, data)
    X <- model.matrix(formula_count, data)
    Y <- model.response(model.frame(formula_count, data))

    pi <- 1 / (1 + exp(-Z %*% fit$coefficients$zero))
    mu <- exp(X %*% fit$coefficients$count)
    theta <- fit$theta
    prob_nb <- (1 - pi) * dnbinom(Y, mu = mu, size = theta)
    prob_zero <- pi * (Y == 0)
    tau <- prob_zero / (prob_zero + prob_nb)
    pred_source <- ifelse(tau > 0.5, "zero_part", "count_part")
    table <- table(pred_source, data$source)
    accuracy <- sum(pred_source == data$source) / length(pred_source)
    return(list(pred_source = pred_source, table = table, accuracy = accuracy))
}
res3 <- predict2(fit3, formula_count = Y ~ X1, formula_zero = Y ~ Z1 + Z2, data)
res$table

res <- predict(fit, formula_count = Y ~ X1, formula_zero = Y ~ Z1 + Z2, data)
res$table
res$accuracy

compare_zinb <- function(beta, gamma, n = 2000) {
    set.seed(42)
    data <- generate_zinb_data(2000, beta, gamma, theta = 2)
    fit1 <- ZINB_Reg(formula_count = Y ~ X1,
                     formula_zero = Y ~ Z1 + Z2,
                     data = data,
                     type = "bWald")
    fit2 <- zeroinfl(Y ~ X1 | Z1 + Z2, data = data, dist = "negbin")
    fit3 <- summary(glmmTMB(formula = Y ~ X1,  ziformula = ~ Z1 + Z2, family = nbinom2, data = data))

    bias1_count <- mean(fit1$coefficients_count - beta)
    bias1_sd_count <- sd(fit1$coefficients_count - beta)
    bias1_zero <- mean(fit1$coefficients_zero - gamma)
    bias1_sd_zero <- sd(fit1$coefficients_zero - gamma)

    bias2_count <- mean(fit2$coefficients$count - beta)
    bias2_sd_count <- sd(fit2$coefficients$count - beta)
    bias2_zero <- mean(fit2$coefficients$zero - gamma)
    bias2_sd_zero <- sd(fit2$coefficients$zero - gamma)

    bias3_count <- mean(fit3$coefficients$cond[, 1] - beta)
    bias3_sd_count <- sd(fit3$coefficients$cond[, 1] - beta)
    bias3_zero <- mean(fit3$coefficients$zi[, 1] - gamma)
    bias3_sd_zero <- sd(fit3$coefficients$zi[, 1] - gamma)

    sd1_count <- mean(fit1$bwald$count$SE)
    sd1_zero <- mean(fit1$bwald$zero$SE)

    sd2_count <- mean(summary(fit2)$coefficients$count[, 2])
    sd2_zero <- mean(summary(fit2)$coefficients$zero[, 2])

    sd3_count <- mean(fit3$coefficients$cond[, 2])
    sd3_zero <- mean(fit3$coefficients$zi[, 2])

    rmse1_count <- sqrt(mean((fit1$coefficients_count - beta)^2))
    rmse1_zero <- sqrt(mean((fit1$coefficients_zero - gamma)^2))

    rmse2_count <- sqrt(mean((fit2$coefficients$count - beta)^2))
    rmse2_zero <- sqrt(mean((fit2$coefficients$zero - gamma)^2))

    rmse3_count <- sqrt(mean((fit3$coefficients$cond[, 1] - beta)^2))
    rmse3_zero <- sqrt(mean((fit3$coefficients$zi[, 1] - gamma)^2))    

    return(list(
        Bias = c(Alg_bWald_count = bias1_count,
                 Alg_bWald_zero = bias1_zero,
                 Alg_zein_count = bias2_count,
                 Alg_zein_zero = bias2_zero,
                 Alg_tmb_count = bias3_count,
                 Alg_tmb_zero = bias3_zero),
        Bias_SD = c(Alg_bWald_count = bias1_sd_count,
                    Alg_bWald_zero = bias1_sd_zero,
                    Alg_zein_count = bias2_sd_count,
                    Alg_zein_zero = bias2_sd_zero,
                    Alg_tmb_count = bias3_sd_count,
                    Alg_tmb_zero = bias3_sd_zero),
        SD = c(Alg_bWald_count = sd1_count,
               Alg_bWald_zero = sd1_zero,
               Alg_zein_count = sd2_count,
               Alg_zein_zero = sd2_zero,
               Alg_tmb_count = sd3_count,
               Alg_tmb_zero = sd3_zero),
        RMSE = c(Alg_bWald_count = rmse1_count,
               Alg_bWald_zero = rmse1_zero,
               Alg_zein_count = rmse2_count,
               Alg_zein_zero = rmse2_zero,
               Alg_tmb_count = rmse3_count,
               Alg_tmb_zero = rmse3_zero)
    ))
}
res4 <- compare_zinb(beta = c(2, -0.5), gamma = c(-1, 1.5, -0.5))

res4 <- compare_zinb(beta = c(2, -0.3, -0.5), gamma = c(-1, 1.5, -0.5))

xtable(as.data.frame(res4), digits = 3)
library(xtable)

indices <- which(res$pred_source == 'zero_part')
zero <- data[indices,]
nb <- data[-indices,]

f1 <- 
    ggplot(data, aes(x = Y)) +
    geom_histogram(aes(y = ..density..), binwidth = diff(range(Y)) / 34, fill = "skyblue", color = "white") +
    labs(x = "Value", y = "Density") +
    theme_minimal()

f2 <- 
    ggplot(nb, aes(x = Y)) +
    geom_histogram(aes(y = ..density..), binwidth = diff(range(Y)) / 34, fill = "skyblue", color = "white") +
    labs(x = "Value", y = "Density") +
    theme_minimal()

grid.arrange(f1, f2, ncol = 2)

# 无随机效应无法对比
library(NBZIMM)
fit5 <- glmm.zinb(
  fixed = Y ~ X1,                   
  random = ~1,
  data = data,                     
  zi_fixed = ~ Z1 + Z2          
)



