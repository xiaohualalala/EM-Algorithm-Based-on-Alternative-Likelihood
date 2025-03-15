library(Rcpp)
library(VGAM)
library(MASS)
library(maxLik)
library(pscl)
library(gamlss)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggrepel)
library(pheatmap)

setwd("/Users/guosa/Desktop/日常/科研/负二项分布/code/NB-EM-Algorithm/Simulation")
source("../Reg/NB_Reg.r", chdir = TRUE)
source("../Reg/ZINB_Reg.r", chdir = TRUE)
source("../Simulation/Dealem.r", chdir = TRUE)
source("../Simulation/Desingle.r", chdir = TRUE)
source("../Simulation/Dtype.r", chdir = TRUE)

sourceCpp("../MLE/C/ZINB_MLE.cpp")

norm <- function(counts) {
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
    return(as.data.frame(counts_norm))
}

# 数据处理
counts <- read.csv("../Simulation/ZINB/counts.csv", row.names = 1)
colnames(counts) <- sub("\\..*", "", colnames(counts))
counts <- cbind(counts[, colnames(counts) == "E3"], counts[, colnames(counts) == "E4"])

# 零百分比
zero_perc <- rowSums(counts == 0) / ncol(counts) * 100
zero_data <- data.frame(zero_perc = zero_perc)
f1 <- ggplot(zero_data, aes(x = zero_perc)) +
    geom_histogram(binwidth = 4, color = "white", fill = "skyblue") + # 修改 binwidth 可调整直方图分辨率
    theme_minimal() +
    labs(
        x = "Zero Percentage (%)",
        y = "Number of Genes"
    )
ggsave("/Users/guosa/Desktop/毕业论文/figures/zero_perc.pdf", f1)


group <- factor(c(rep(1, 81), rep(2, 190)))

# Detecting the differentially expressed genes
results <- DEsingle(counts = counts, group = group)
results2 <- DEalem(counts = counts, group = group)
results2_valid <- results2[is.na(results2$Remark), ]

table(results$Remark)
table(results2$Remark)

# Dividing the differentially expressed genes into 3 categories
results.classified <- DEtype(results = results, threshold = 0.05)
results2.classified <- DEtype(results = results2_valid, threshold = 0.05)

table(results.classified$Type)
table(results2.classified$Type)

DEg <- subset(results.classified, Type == "DEg")
DEs <- subset(results.classified, Type == "DEs")
DEa <- subset(results.classified, Type == "DEa")

DEg2 <- subset(results2.classified, Type == "DEg")
DEs2 <- subset(results2.classified, Type == "DEs")
DEa2 <- subset(results2.classified, Type == "DEa")

counts_norm <- norm(counts)

DEs["ICAM5", ]
DEs2["ICAM5", ]

counts_I <- counts_norm %>%
    filter(rownames(counts_norm) == "ICAM5")

counts_I <- data.frame(
    group = c(rep("E3", 81), rep("E4", 190)), # 分组信息
    value = as.numeric(counts_I)
)

f2 <- ggplot(counts_I, aes(x = value, fill = group)) +
    geom_histogram(
        binwidth = diff(range(counts_I$value)) / 35,
        position = "identity", # 允许不同组的直方图重叠
        alpha = 0.7 # 设置透明度
    ) +
    scale_fill_manual(values = c("E3" = "salmon", "E4" = "skyblue")) + # 自定义颜色
    theme_minimal() +
    labs(
        title = "DEs gene ICAM5",
        x = "Normalized read counts",
        y = "Number of cells",
        fill = "group"
    )


DEa["EIF3M", ]
DEa2["EIF3M", ]

counts_E <- counts_norm %>%
    filter(rownames(counts_norm) == "EIF3M")

counts_E <- data.frame(
    group = c(rep("E3", 81), rep("E4", 190)), # 分组信息
    value = as.numeric(counts_E)
)

f3 <- ggplot(counts_E, aes(x = value, fill = group)) +
    geom_histogram(
        binwidth = diff(range(counts_E$value)) / 35,
        position = "identity", # 允许不同组的直方图重叠
        alpha = 0.7 # 设置透明度
    ) +
    scale_fill_manual(values = c("E3" = "salmon", "E4" = "skyblue")) + # 自定义颜色
    theme_minimal() +
    labs(
        title = "DEa gene HSPA5",
        x = "Normalized read counts",
        y = "Number of cells",
        fill = "group"
    )

DEg["TMEM14B", ]
DEg2["TMEM14B", ]

counts_T <- counts_norm %>%
    filter(rownames(counts_norm) == "TMEM14B")

counts_T <- data.frame(
    group = c(rep("E3", 81), rep("E4", 190)), # 分组信息
    value = as.numeric(counts_T)
)


f4 <- ggplot(counts_T, aes(x = value, fill = group)) +
    geom_histogram(
        binwidth = diff(range(counts_T$value)) / 35,
        position = "identity", # 允许不同组的直方图重叠
        alpha = 0.7 # 设置透明度
    ) +
    scale_fill_manual(values = c("E3" = "salmon", "E4" = "skyblue")) + # 自定义颜色
    theme_minimal() +
    labs(
        title = "DEg gene TMEM14B",
        x = "Normalized read counts",
        y = "Number of cells",
        fill = "group"
    )

f5 <- grid.arrange(f2, f3, f4, nrow = 3)
ggsave("/Users/guosa/Desktop/毕业论文/figures/DeGene.pdf", f5, width = 5, height = 12.5, dpi = 300)

# 火山图

results2.classified$Gene <- row.names(results2.classified)
results2.classified$log2FoldChange <- log2(results2.classified$norm_foldChange)

# 读取火山图数据文件 选取3000个基因
set.seed(432)
data <- subset(results2.classified, pvalue.adj.FDR > 0 & abs(log2FoldChange) < 10)
index <- sample(1:15602, 8000, replace = FALSE)
data <- data[index, ]

# 设置p_value和logFC的阈值
cut_off_fdr <- 0.05 # 统计显著性
cut_off_logFC <- 1 # 差异倍数值

# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
data$sig <- ifelse(data$pvalue.adj.FDR < cut_off_fdr & abs(data$log2FoldChange) >= cut_off_logFC,
    ifelse(data$log2FoldChange > cut_off_logFC, "Up", "Down"),
    "Not Sig"
)

# 根据PValue小于多少和log[2]FC的绝对值大于多少筛选出合适的点
data$label <- ifelse(data$pvalue.adj.FDR < 0.0000001 & abs(data$log2FoldChange) > 5,
    as.character(data$Gene), ""
)

table(data$sig)

f6 <- ggplot(data, aes(data$log2FoldChange, -log10(data$pvalue.adj.FDR))) + # 定义横纵坐标
    geom_point(alpha = 0.8, size = 1.5, aes(color = sig)) + # 绘制散点图，分组依据是数据框的sig列
    labs(
        x = "log[2](Fold Change)",
        y = "-log[10](FDR)"
    ) +
    scale_color_manual(values = c("blue", "#d2dae2", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) + # 在图上添加虚线
    geom_vline(xintercept = c(-1, 1), linetype = 2) + # 在图上添加虚线
    geom_text_repel(
        aes(
            x = data$log2FoldChange, # geom_text_repel 标记函数
            y = -log10(data$pvalue.adj.FDR),
            label = label
        ),
        max.overlaps = 50, # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
        size = 3, # 字体大小
        box.padding = unit(0.5, "lines"), # 标记的边距
        point.padding = unit(0.1, "lines"),
        segment.color = "darkgrey", # 标记线条的颜色
        show.legend = FALSE
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank()
    )


ggsave("/Users/guosa/Desktop/毕业论文/figures/volcano.pdf", f6)

up_gene <- subset(data, Gene == 'ZNF25')

log10(up_gene$pvalue.adj.FDR)
up_gene$log2FoldChange
2^up_gene$log2FoldChange


# 聚类图
# 每种类型基因的表达水平

DEg_names <- row.names(DEg2)
counts_norm_DEg <- subset(counts_norm, rownames(counts_norm) %in% DEg_names)

DEs_names <- row.names(DEs2)
counts_norm_DEs <- subset(counts_norm, rownames(counts_norm) %in% DEs_names)

DEa_names <- row.names(DEa2)
counts_norm_DEa <- subset(counts_norm, rownames(counts_norm) %in% DEa_names)

heatmatrix <- rbind(rbind(counts_norm_DEg[1:400, ], counts_norm_DEs[1:400, ]), counts_norm_DEa[1:400, ])

heatmatrix <- log(heatmatrix + 1)

# 行注释信息
annotation_row <- data.frame(
    GeneClass = rep(c("DEg", "DEs", "DEa"), each = 400)
)
rownames(annotation_row) <- rownames(heatmatrix)

# 列注释信息
annotation_col <- data.frame(
    CellType = c(rep("E3", each = 81), rep("E4", each = 190))
)
rownames(annotation_col) <- colnames(heatmatrix)

annotation_colors <- list(
    GeneClass = c(
        "DEg" = "#1F77B4", # 明亮蓝色
        "DEs" = "#2CA02C", # 鲜亮绿色
        "DEa" = "#FFCC00"  # 亮黄色
    )
)

# 第一步：先用 silent = TRUE 运行一次 pheatmap 得到聚类结果
res <- pheatmap(heatmatrix,
    gaps_col = 81,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    clustering_distance_rows = "euclidean", # 使用欧氏距离
    clustering_method = "complete", # 使用完全链接法
    cluster_rows = TRUE, # 启用行聚类
    cluster_cols = FALSE, # 禁用列聚类（可选）
    silent = TRUE, # 静默模式，不直接绘图
    border = FALSE,
    annotation_colors = annotation_colors
)

# 第二步：根据聚类后的行顺序生成标签
# 获取聚类结果返回的行顺序 (res$tree_row$order 是一个整数向量)
ordered_row_names <- rownames(heatmatrix)[res$tree_row$order]

selected_genes <- ordered_row_names[seq(1, 1200, 60)]

selected_cells <- colnames(heatmatrix)[c(seq(1, 81, 9), seq(82, 271, 9))]

# 第三步：重新绘制热图，并使用自定义的 labels_row 参数
f7 <- pheatmap(heatmatrix,
    gaps_col = 81,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    clustering_distance_rows = "euclidean", # 使用欧氏距离
    clustering_method = "complete", # 使用完全链接法
    cluster_rows = TRUE, # 启用行聚类
    cluster_cols = FALSE, # 禁用列聚类（可选）
    labels_row = ifelse(rownames(heatmatrix) %in% selected_genes,
        rownames(heatmatrix), ""
    ), # 使用定制的行标签
    labels_col = ifelse(colnames(heatmatrix) %in% selected_cells,
        colnames(heatmatrix), ""
    ),
    border = FALSE,
    annotation_colors = annotation_colors,
    fontsize = 7,
    fontsize_row = 7,
    fontsize_col = 7,
    legend = TRUE,
    annotation_legend = TRUE,
    angle_col = 45 # 列注释倾斜 45 度
)

ggsave("/Users/guosa/Desktop/毕业论文/figures/heatmap.pdf", f7, width = 10, height = 8, dpi = 300)
