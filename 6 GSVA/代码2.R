rm(list = ls())
setwd("D:\\科研\\王晓 宫颈癌 铜死亡分型 结果+组图\\6 GSVA")
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

# 读取表达数据并预处理 ------------------------------------------------
load("merge.RDATA")
data <- as.matrix(outTab)
data <- avereps(data)

# 读取基因集 --------------------------------------------------------
geneSets <- getGmt("c2.cp.reactome.v7.5.1.symbols.gmt")

# GSVA 分析 --------------------------------------------------------
param <- gsvaParam(exprData = data, geneSets = geneSets)
gsva_es <- gsva(param)
write.table(gsva_es, file = "gsvaOut1.txt", sep = "\t", quote = F, col.names = NA)

# 读取分型信息并合并 ------------------------------------------------
cluster <- read.table("Cluster.txt", header = T, sep = "\t", row.names = 1)
cluster$cluster <- factor(cluster$cluster)  # 分型为 A/B

# 确保样本名匹配并提取数据
sameSample <- intersect(colnames(gsva_es), rownames(cluster))
gsva_es_subset <- gsva_es[, sameSample, drop = FALSE]
cluster_subset <- cluster[sameSample, , drop = FALSE]

# 构建分析数据框（添加 Project 列）-----------------------------------
analysis_df <- data.frame(
  t(gsva_es_subset),  # 转置为 [样本 x 通路]
  cluster = cluster_subset$cluster,
  Project = gsub("(.*?)\\_.*", "\\1", rownames(cluster_subset))  # 从样本名提取 Project
)

# 定义比较组（A vs B）------------------------------------------------
allType <- levels(analysis_df$cluster)
comp <- combn(allType, 2)

# 差异分析循环 ----------------------------------------------------------
adj.P.Val.Filter <- 0.05
bioCol <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")

for (i in 1:ncol(comp)) {
  # 样本分组
  treat <- analysis_df[analysis_df$cluster == comp[2,i], ]
  con <- analysis_df[analysis_df$cluster == comp[1,i], ]
  combined_data <- rbind(con, treat)
  
  # 提取纯数值矩阵
  expr_matrix <- as.matrix(combined_data[, !colnames(combined_data) %in% c("cluster", "Project")])
  
  # 构建设计矩阵
  design <- model.matrix(~0 + combined_data$cluster)
  colnames(design) <- levels(combined_data$cluster)
  
  # 差异分析
  fit <- lmFit(t(expr_matrix), design)
  contrast <- paste0(comp[2,i], "-", comp[1,i])
  cont_matrix <- makeContrasts(contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)
  
  # 输出所有差异通路
  allDiff <- topTable(fit2, adjust = "fdr", number = Inf)
  allDiffOut <- rbind(id = colnames(allDiff), allDiff)
  write.table(allDiffOut, file = paste0(contrast, ".all.txt"), sep = "\t", quote = F, col.names = F)
  
  # 输出显著差异通路（关键修复：确保 diffSig 在条件判断内定义）
  if (exists("allDiff") && nrow(allDiff) > 0) {
    diffSig <- allDiff[with(allDiff, (abs(logFC) > 0.1 & adj.P.Val < adj.P.Val.Filter)), ]
    if (nrow(diffSig) > 0) {
      diffSigOut <- rbind(id = colnames(diffSig), diffSig)
      write.table(diffSigOut, file = paste0(contrast, ".diff.txt"), sep = "\t", quote = F, col.names = F)
      
      # 热图绘制（仅在 diffSig 存在时执行）
      # 准备注释数据（包含 cluster 和 Project）
      ann <- data.frame(
        cluster = combined_data$cluster,
        Project = combined_data$Project
      )
      rownames(ann) <- rownames(combined_data)
      
      # 设置颜色
      CluCol <- bioCol[1:length(levels(analysis_df$cluster))]
      names(CluCol) <- levels(analysis_df$cluster)
      CluCol <- c("A" = "#0066FF", "B" = "#FF9900")
      ProjCol <- c("Project1" = "#ff9389", "Project2" = "#7CC767")
      names(ProjCol) <- unique(analysis_df$Project)
      
      ann_colors <- list(
        cluster = CluCol[comp[,i]],
        Project = ProjCol
      )
      
      # 提取热图数据
      termNum <- 20
      diffTermName <- rownames(diffSig)
      diffLength <- length(diffTermName)
      if (diffLength < termNum) termNum <- diffLength
      hmGene <- diffTermName[1:termNum]
      hmExp <- expr_matrix[, hmGene, drop = FALSE]
      # 修改热图通路名，加上 q 值
      annotated_names <- paste0(hmGene, "\n(q=", signif(allDiff[hmGene, "adj.P.Val"], 2), ")")
      colnames(hmExp) <- annotated_names
      
      # 绘制热图
      pdf(file = paste0(contrast, "reactome.heatmap.pdf"), height = 6, width = 10)
      pheatmap(
        t(hmExp),
        annotation_col = ann,
        annotation_colors = ann_colors,
        color = colorRampPalette(c(rep("blue",2), "black", rep("yellow",2)))(50),
        cluster_cols = FALSE,
        show_colnames = FALSE,
        gaps_col = cumsum(table(combined_data$cluster)),
        scale = "row",
        fontsize = 10,
        fontsize_row = 7,
        fontsize_col = 10
      )
      dev.off()
    }
  }
}

