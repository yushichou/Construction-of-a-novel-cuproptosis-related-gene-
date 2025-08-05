rm(list = ls())
setwd("D:\\科研\\王晓 宫颈癌 铜死亡分型 结果+组图\\21 score_基因表达\\new")

# 加载包
library(limma)
library(ggpubr)

# 文件路径
scoreFile <- "score.group.txt"
load("merge.RDATA")
data <- outTab
colnames(data) <- gsub("(.*?)\\_(.*?)", "\\2", colnames(data))
data2 <- data

# 读取打分分组文件
score <- read.table(scoreFile, header = T, sep = "\t", check.names = F, row.names = 1)

# 定义函数
plot_gene_boxplot_with_qval <- function(gene, expr_data, score_data, output_dir = ".", palette_colors = c("#00AFBB", "#E7B800")) {
  # 提取目标基因表达
  if (!gene %in% rownames(expr_data)) {
    warning(paste("Gene", gene, "not found in expression matrix. Skipped."))
    return(NULL)
  }
  
  gene_expr <- rbind(expr_data, gene = expr_data[gene, ])
  exp <- t(gene_expr[c("gene", gene), ])
  exp <- avereps(exp)
  
  # 合并表达数据和打分
  sameSample <- intersect(rownames(exp), rownames(score_data))
  exp <- exp[sameSample, ]
  score_sub <- score_data[sameSample, ]
  df <- cbind(as.data.frame(exp), as.data.frame(score_sub))
  df$group <- factor(df$group, levels = c("Low", "High"))
  
  # Wilcoxon 检验 + q 值
  wilcox_res <- wilcox.test(gene ~ group, data = df)
  p_val <- wilcox_res$p.value
  q_val <- p.adjust(p_val, method = "fdr")
  qval_label <- paste0("q = ", formatC(q_val, format = "e", digits = 2))
  
  # 绘图
  boxplot <- ggboxplot(df, x = "group", y = "gene", fill = "group",
                       xlab = "Score Group",
                       ylab = paste(gene, "expression"),
                       legend.title = "Score",
                       palette = palette_colors) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Low", "High"))) +
    annotate("text", x = 1.5, y = max(df$gene) * 1.05, label = qval_label, size = 5)
  
  # 输出文件
  pdf(file = file.path(output_dir, paste0(gene, "_qval_boxplot.pdf")), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}

# 指定基因
genes_to_plot <- c("CD274", "CTLA4", "LAG3", "PDCD1", "TIGIT")

# 批量绘图
for (gene in genes_to_plot) {
  plot_gene_boxplot_with_qval(gene = gene, expr_data = data2, score_data = score)
}

