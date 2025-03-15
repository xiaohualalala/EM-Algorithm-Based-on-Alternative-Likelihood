library(Rcpp)
library(ggplot2)
require(pscl)
library(xtable)
library(gridExtra)

setwd("/Users/guosa/Desktop/日常/科研/负二项分布/code/NB-EM-Algorithm/Simulation")

sourceCpp("../Simulation/Samples/Generate_Samples.cpp")
sourceCpp("../Simulation/NB/SimulateNB.cpp")
sourceCpp("../Simulation/MixNB/SimulateMixNB.cpp")
sourceCpp("../Simulation/ZINB/SimulateZINB.cpp")
sourceCpp("../MLE/C/ZINB_MLE.cpp")
sourceCpp("../MLE/C/NB_MLE.cpp")
source('../Simulation/MixNB/VisualMixNB.r')

# 设置随机数种子
set.seed(123)

# NB
# 模拟参数
prob_values <- seq(0.2, 0.8, 0.2)
size_values <- seq(1, 9, 2)
res_NB <- SimulateNBs(prob_values, size_values, 200, 2000)

res_NB_latex <- xtable(res_NB)
# 迭代次数与参数的热力图
p1 <- 
ggplot(res_NB, aes(x = as.factor(size), y = as.factor(prob))) +
    geom_tile(aes(fill = avg_t)) +
    scale_fill_gradient(low = "#ffeda0", high = "#f03b20") +
    labs(x = "size", y = "prob", fill = "Iterations") +
    theme_minimal()

ggsave("/Users/guosa/Desktop/毕业论文/figures/hotplot.pdf", p1)

# 不同初始点对收敛次数的影响
source('../MLE/R/NB_MLE.r', chdir = TRUE)

# 生成模拟数据
set.seed(123)
X <- rnbinom(2000, size = 5, prob = 0.5)

# random
random_conver <- logical(200)
random_time <- numeric(200)
random_iter <- numeric(200)

for (i in 1:length(conver)) {
    init_values <- c(runif(1, 0.1, 50), runif(1, 0.01, 0.99))
    result <- NB_MLE(X, init_values = init_values, verbose = FALSE)
    random_conver[i] <- result$converged
    random_time[i] <- result$duration
    random_iter[i] <- result$iter
}

random_rate <- sum(random_conver) / length(random_conver)
random_cpu <- mean(random_time)
sd(random_time)
mean(random_iter)
sd(random_iter)

# extreme 

extreme_conver <- logical(4)
extreme_time <- numeric(4)
extreme_iter <- numeric(4)

init_values <- list(c(0.1,0.001),c(0.1,0.999),c(1000,0.001),c(1000,0.999))

for (i in 1:length(init_values)) {
    result <- NB_MLE(X, init_values = unlist(init_values[i]), verbose = FALSE)
    extreme_conver[i] <- result$converged
    extreme_time[i] <- result$duration
    extreme_iter[i] <- result$iter
}
extreme_rate <- sum(extreme_conver) / length(extreme_conver)
extreme_cpu <- mean(extreme_time)
sd(extreme_time)
mean(extreme_iter)
sd(extreme_iter)

# moment 
result <- NB_MLE(X, verbose = FALSE)
sum(result$converged)
result$duration
result$iter

library(plotly)

# 定义初始值的范围
prob_values <- seq(0.1, 0.9, length.out = 40)  # 取10个概率值
size_values <- seq(1, 20, length.out = 40)       # 取10个size值

init_values <- c(runif(1, 0.1, 0.9), runif(1, 1, 20))
# 记录结果
results <- expand.grid(Prob = prob_values, Size = size_values)
results$Iterations <- NA

for (i in 1:nrow(results)) {
    init_values <- c(results$Size[i], results$Prob[i])
    result <- NB_MLE(X, init_values = init_values, verbose = FALSE)
    results$Iterations[i] <- result$iter
}

# 确保数据格式正确
Iterations <- matrix(results$Iterations, nrow = length(prob_values), ncol = length(size_values))

# 创建三维曲面图
p2 <- 
plot_ly(
    x = ~prob_values,
    y = ~size_values,
    z = ~Iterations,
    type = 'surface'
) %>%
layout(
    scene = list(
        xaxis = list(title = 'Prob'),
        yaxis = list(title = 'Size'),
        zaxis = list(title = 'Iterations')
    )
)
library(htmlwidgets)
saveWidget(p2, "/Users/guosa/Desktop/毕业论文/figures/初始值与迭代次数.html")

## MixNB
source("../Simulation/MixNB/VisualMixNB.r")

# 双样本
# 模拟参数
set.seed(123)
weights <- c(0.4, 0.6)
probs <- c(0.125, 0.25)
sizes <- c(5, 30)

# 生成数据并进行分组
res2 <- SimulateMixNB(probs, sizes, weights, 2000)

# 对比图
p3 <- plotHistDensity(res2$observations, res2$estimation, res2$accuracy)
ggsave("/Users/guosa/Desktop/毕业论文/figures/mix2.pdf", p3)

set.seed(234)
weights <- c(0.4, 0.6)
means <- c(50, 100)
vars <- c(400, 400)
probs <- means / vars
sizes <- means^2 / (vars - means)

res2 <- SimulateMixNB(probs, sizes, weights, 2000)

p3 <- plotHistDensity(res2$observations, res2$estimation, res2$accuracy)

ggsave("/Users/guosa/Desktop/yudabian0221/figures/mix2.pdf", p3)

source('../Simulation/MixNB/VisualMixNB.r')

# 运行100次，观察准确率与迭代次数
data <- SimulateMixNBs(probs, sizes, weights, 100, 2000)

df_accuracy <- data.frame(Value = data$accuracy, Type = "Accuracy")
mean(df_accuracy$Value)
sd(df_accuracy$Value)

df_iterations <- data.frame(Value = data$iterations, Type = "Iterations")
mean(df_iterations$Value)
sd(df_iterations$Value)

K <- length(weights)
df_size <- data$eSize
size_mean <- sapply(1:K, function(i) mean(sapply(df_size, `[`, i)))
size_sd <- sapply(1:K, function(i) sd(sapply(df_size, `[`, i)))

df_prob <- data$eProb
prob_mean <- sapply(1:K, function(i) mean(sapply(df_prob, `[`, i)))
prob_sd <- sapply(1:K, function(i) sd(sapply(df_prob, `[`, i)))

df_weight <- data$eWeight
weight_mean <- sapply(1:K, function(i) mean(sapply(df_weight, `[`, i)))
weight_sd <- sapply(1:K, function(i) sd(sapply(df_weight, `[`, i)))

# accuracy箱线图
p5 <- 
ggplot(df_accuracy, aes(x = Type, y = Value)) +
    geom_boxplot() +
    theme_minimal() +
    ylim(0.965, 0.975) + 
    annotate("text", x = Inf, y = Inf, label = paste("Accuracy:", round(mean(df_accuracy$Value),4),"\u00B1", round(sd(df_accuracy$Value),4)), hjust = 1.1, vjust = 2, size = 8, color = "blue") + 
    labs(x = "", y = "Accuracy")

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix4.pdf", p5)

# iterations箱线图
p6 <- 
ggplot(df_iterations, aes(x = Type, y = Value)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "", y = "Iterations")

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix5.pdf", p6)

# 模拟参数

weights <- c(0.2, 0.3, 0.5)
probs <- c(0.1, 0.3, 0.2)
sizes <- c(10, 125, 50)

# 生成数据并进行分组
res3 <- SimulateMixNB(probs, sizes, weights, 2000)

# 对比图
p4 <- plotHistDensity(res3$observations, res3$estimation, res3$accuracy)

ggsave("/Users/guosa/Desktop/毕业论文/figures/mix3.pdf", p4)

set.seed(234)
weights <- c(0.2, 0.3, 0.5)
means <- c(100, 300, 200)
vars <- c(1000, 1000, 1000)
probs <- means / vars
sizes <- means^2 / (vars - means)

res3 <- SimulateMixNB(probs, sizes, weights, 2000)

p4 <- plotHistDensity(res3$observations, res3$estimation, res3$accuracy)

ggsave("/Users/guosa/Desktop/yudabian0221/figures/mix3.pdf", p4)

# 运行100次，观察准确率与迭代次数
data <- SimulateMixNBs(probs, sizes, weights, 100, 2000)

df_accuracy <- data.frame(Value = data$accuracy, Type = "Accuracy")
mean(df_accuracy$Value)
sd(df_accuracy$Value)
df_iterations <- data.frame(Value = data$iterations, Type = "Iterations")
mean(df_iterations$Value)
sd(df_iterations$Value)

K <- length(weights)
df_size <- data$eSize
size_mean <- sapply(1:K, function(i) mean(sapply(df_size, `[`, i)))
size_sd <- sapply(1:K, function(i) sd(sapply(df_size, `[`, i)))

df_prob <- data$eProb
prob_mean <- sapply(1:K, function(i) mean(sapply(df_prob, `[`, i)))
prob_sd <- sapply(1:K, function(i) sd(sapply(df_prob, `[`, i)))

df_weight <- data$eWeight
weight_mean <- sapply(1:K, function(i) mean(sapply(df_weight, `[`, i)))
weight_sd <- sapply(1:K, function(i) sd(sapply(df_weight, `[`, i)))

p5 <- grid.arrange(p3, p4, nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/mix6.pdf", p5, width = 8, height = 7.5, dpi = 300)

ggsave("/Users/guosa/Desktop/yudabian0221/figures/mix6.pdf", p5)

# confusion matrix
conf_matrix <- confusionMatrix(as.factor(res2$predicted_labels), as.factor(res2$true_labels))
conf_matrix

# roc
library(pROC)
library(caret)

roc_curve2 <- roc(res2$true_labels, res2$postProbs[,1], plot = TRUE, print.auc = TRUE)

roc_curve3_1 <- roc(res3$true_labels == 1, res3$postProbs[,1], plot = TRUE, print.auc = TRUE)

roc_curve3_2 <- roc(res3$true_labels == 2, res3$postProbs[,2], plot = TRUE, print.auc = TRUE)

roc_curve3_3 <- roc(res3$true_labels == 3, res3$postProbs[,3], plot = TRUE, print.auc = TRUE)

pdf("/Users/guosa/Desktop/毕业论文/figures/ROC_Curve.pdf")
plot(roc_curve2, col = "coral", main = "ROC Curves", print.auc = TRUE, auc.polygon = FALSE)
plot(roc_curve3_3, col = "navy", add = TRUE, print.auc = TRUE, auc.polygon = FALSE, print.auc.y = 0.4)
legend("bottomright", legend = c("K=2", "K=3"), col = c("coral", "navy"), lwd = 2)
dev.off()


roc_list <- list()
auc_list <- list()



# 遍历每个类别
for (i in 1:3) {
  # 提取当前类别的真实标签（二值化）
  true_labels <- as.numeric(res3$true_labels == levels(factor(res3$true_labels)[i]))
  
  # 计算ROC曲线和AUC
  roc_obj <- roc(true_labels, res3$postProbs[, i])
  auc_value <- auc(roc_obj)
  
  # 存储结果
  roc_list[[i]] <- roc_obj
  auc_list[[i]] <- auc_value
}
# 创建绘图数据框
plot_data <- data.frame()
for (i in 1:3) {
  temp_df <- data.frame(
    FPR = 1 - roc_list[[i]]$specificities,
    TPR = roc_list[[i]]$sensitivities,
    Class = paste("Class", i-1, " (AUC=", round(auc_list[[i]], 2), ")", sep="")
  )
  plot_data <- rbind(plot_data, temp_df)
}

# 绘制曲线
ggplot(plot_data, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "False Positive Rate", y = "True Positive Rate", 
       title = "ROC Curves for 3-Class Classification") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red", "green"))


install.packages("pROC") 
install.packages("multiROC")
library(multiROC)
library(pROC)



## ZINB
# 模拟参数
prob_values <- seq(0.2, 0.6, 0.2)
size_values <- seq(6, 10, 2)
zero_values <- seq(0.2, 0.6, 0.2)

res_ZINB <- SimulateZINBs(prob_values, size_values, zero_values, 20, 2000)
res_ZINB_latex <- xtable(res_ZINB)

# 分类效果模拟
set.seed(123)
prob <- 0.1
size <- 2
zero <- 0.2

ZINB <- SamplesZINB(prob, size, zero, 1, 2000)
X <- unlist(ZINB$samples)
true_labels <- unlist(ZINB$labels)

res_ZINB <- ZINB_MLE(X)
eZero <- res_ZINB$Parameter_Estimates$eZero
eSize <- res_ZINB$Parameter_Estimates$eSize
eProb <- res_ZINB$Parameter_Estimates$eProb
prob_nb <- (1 - eZero) * dnbinom(X, size = eSize, prob = eProb)
prob_zero <- eZero * (X == 0)
tau <- prob_zero / (prob_zero + prob_nb)
predict_labels <- ifelse(tau > 0.5, 0, 1)

eZero / ((1 - eZero) * dnbinom(0, size = eSize, prob = eProb) + eZero)
sum(X==0) dnbinom(0, size = eSize, prob = eProb)*2000

accuracy <- mean(true_labels == predict_labels)

NB_part <- X[predict_labels == 1]
Zero_part <- X[predict_labels == 0]
Source <- ifelse(tau > 0.5, "Negative Binomial", "Zero-Inflated")

# 计算估计的负二项分布的概率密度
x_values <- seq(0, max(NB_part), by = 1)
nb_density <- dnbinom(x_values, size = eSize, prob = eProb)

# 创建数据框用于绘制概率密度曲线
density_data <- data.frame(x_values, nb_density)

data <- data.frame(Value = X)

# 绘制直方图
p7 <- 
    ggplot(data, aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = diff(range(X)) / 70, fill = "skyblue", color = "white") +
    labs(x = "Value",
    y = "Density") +
    theme_minimal() 

ggsave("/Users/guosa/Desktop/毕业论文/figures/zinb-hist.pdf", p7)

  +
  annotate("text", x = Inf, y = Inf, label = paste("Accuracy:", accuracy), hjust = 1.1, vjust = 2, size = 5, color = "blue")


# 添加概率密度曲线
p8 <- 
p7 + geom_line(data = density_data, aes(x = x_values, y = nb_density), color = "red", size = 1)

ggsave("/Users/guosa/Desktop/毕业论文/figures/zinb1.pdf", p8)

# Rate 
set.seed(123)

probs <- seq(0.1, 0.9, 0.2)
sizes <- seq(1, 20, 2)
rates <- SimulateRates(probs,sizes,20,2000)

p9 <- 
ggplot(rates, aes(x = size, y = avg_rate, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_rate - sd_rate, ymax = avg_rate + sd_rate), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Size", 
       y = "Convergence Rate", 
       color = "Prob", 
       title = "Convergence Rate vs Size") +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate1.pdf", p9)

probs <- seq(0.1, 0.9, 0.05)
sizes <- seq(1, 10, 2)
rates <- SimulateRates(probs,sizes,20,2000)
 
p10 <- 
ggplot(rates, aes(x = prob, y = avg_rate, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  ylim(0.42, 0.98) +
  geom_errorbar(aes(ymin = avg_rate - sd_rate, ymax = avg_rate + sd_rate), width = 0.01, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Prob", y = "Convergence Rate", color = "Size", title = "Convergence Rate vs Prob") +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate2.pdf", p10)

p11 <- grid.arrange(p9, p10, ncol = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/rates.pdf", p11, width = 16, height = 6, dpi = 300)

prob <- 0.6
size <- 1
sample_sizes <- seq(500,100000,500)
rates <- data.frame(n = sample_sizes, rate = numeric(200))
for (i in 1:200){
    rates$rate[i] <- SimulateRate(prob,size,1,sample_sizes[i])$conv_rates
}

p12 <- 
ggplot(rates, aes(x = n, y = rate)) +
  geom_point(color='red', size=0.8) +
  labs(x = "n", y = "Convergence Rate", title= "Convergence Rate vs Sample Size") +
  theme_bw()+
  ylim(0.5, 0.65) +
  annotate("text", x = max(rates$n), y = 0.65, label = "Size=1.0 \n Prob=0.6", hjust=1, vjust=1, color = "blue", size = 5)

ggsave("/Users/guosa/Desktop/毕业论文/figures/rate3.pdf", p12, width = 8, height = 6, dpi = 300)

library(Rcpp)
sourceCpp("../Simulation/NB/SimulateNB.cpp")

set.seed(123)

probs <- seq(0.1, 0.9, 0.2)
sizes <- seq(1, 20, 2)
errors <- SimulateErrors(probs,sizes,20,2000)

p13 <- 
ggplot(errors, aes(x = size, y = avg_error_1, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_1 - sd_error_1, ymax = avg_error_1 + sd_error_1), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.05) + 
  labs(x = "Size", 
       y = bquote("Standard Error (" ~ phi[1] ~ ")"), 
       color = "Prob", 
       title = expression("Standard Error (" ~ phi[1] ~ ") vs Size" )) +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error1.pdf", p13)

p14 <- 
ggplot(errors, aes(x = size, y = avg_error_2, color = as.factor(prob))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.1, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 1) + 
  labs(x = "Size", 
       y = bquote("Standard Error (" ~ phi[2] ~ ")"), 
       color = "Prob",
       title = expression("Standard Error (" ~ phi[2] ~ ") vs Size" )) +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error2.pdf", p14)

probs <- seq(0.1, 0.9, 0.05)
sizes <- seq(1, 10, 2)
errors <- SimulateErrors(probs,sizes,20,2000)

p15 <- 
ggplot(errors, aes(x = prob, y = avg_error_1, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_1 - sd_error_1, ymax = avg_error_1 + sd_error_1), width = 0.01,  linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.05) + 
  labs(x = "Prob", 
       y = bquote("Standard Error (" ~ phi[1] ~ ")"), 
       color = "Size",
       title = expression("Standard Error (" ~ phi[1] ~ ") vs Prob" )) +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error3.pdf", p15)

p16 <- 
ggplot(errors, aes(x = prob, y = avg_error_2, color = as.factor(size))) +
  geom_line(size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_error_2 - sd_error_2, ymax = avg_error_2 + sd_error_2), width = 0.01, linetype = "twodash") +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.5) + 
  labs(x = "Prob", 
       y = bquote("Standard Error (" ~ phi[2] ~ ")"), 
       color = "Size",
       title = expression("Standard Error (" ~ phi[2] ~ ") vs Prob" )) +
  theme_bw()

ggsave("/Users/guosa/Desktop/毕业论文/figures/error4.pdf", p16)

p17 <- grid.arrange(p13,p15,nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/error_1.pdf", p17)

p18 <- grid.arrange(p14,p16,nrow = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/error_2.pdf", p18)

p19 <- grid.arrange(p13,p15,p14,p16,nrow = 2, ncol = 2)
ggsave("/Users/guosa/Desktop/毕业论文/figures/error_3.pdf", p19, width = 16, height = 12, dpi=400)

source("../Simulation/scRNA.r", chdir = TRUE)

