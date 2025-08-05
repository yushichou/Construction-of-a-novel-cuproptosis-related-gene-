# ==== 设置工作环境 ====
rm(list = ls())
setwd("D:\\科研\\王晓 宫颈癌 铜死亡分型 结果+组图\\16 PCA评分\\内部验证3")

# ==== 加载依赖包 ====
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)
library(caret)
library(maxstat)
library(stringr)

# ==== 读取表达矩阵和生存数据 ====
expFile = "uniSigGeneExp.txt"
cliFile = "clinicaldata.txt"
exp = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(exp) <- sub("^TCGA_", "", colnames(exp))
cli = read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime = cli$futime / 365
cli$M_Sstage <- as.factor(cli$M_Sstage)
cli$N_stage <- as.factor(cli$N_stage)
cli$T_stage <- as.factor(cli$T_stage)
cli$New_tumor_event <- factor(cli$New_tumor_event, 
                              levels = c("NO", "YES"), 
                              labels = c("No", "Yes"))
cli$Grade <- as.factor(cli$Grade)
cli$Stage <- as.factor(cli$Stage)

# ==== PCA评分打分 ====
exp = t(exp)
pca <- prcomp(exp, scale=TRUE)
value = predict(pca)
score_raw = value[,1] + value[,2]
score_df <- data.frame(score = as.numeric(score_raw))
rownames(score_df) <- rownames(exp)
write.table(score_df, file="score.txt", sep="\t", quote=F, col.names=F)

# ==== 合并数据 ====
sameSample = intersect(rownames(score_df), rownames(cli))
data_all = cbind(cli[sameSample, ], score_df[sameSample, , drop=FALSE])

# ==== 使用maxstat确定cutoff ====
maxstat_obj <- maxstat.test(Surv(futime, fustat) ~ score, data=data_all, 
                            smethod="LogRank", pmethod="exactGauss")
cutoff <- maxstat_obj$estimate
cat("Optimal cutoff from maxstat:", cutoff, "\n")

# ==== Bootstrap验证函数 ====
bootstrap_validation <- function(data_all, n_boot = 100, times = c(1, 3, 5)) {
  auc_matrix <- matrix(NA, nrow = n_boot, ncol = length(times))
  colnames(auc_matrix) <- paste0(times, "-year")
  pval_vec <- numeric(n_boot)
  roc_list <- vector("list", length = n_boot)
  valid_count <- 0
  
  for (i in 1:n_boot) {
    set.seed(i)
    idx <- sample(1:nrow(data_all), size = floor(0.6 * nrow(data_all)))
    train_data <- data_all[idx, ]
    test_data <- data_all[-idx, ]
    
    if (nrow(test_data) < 20 || nrow(train_data) < 30) next
    
    train_data$group <- ifelse(train_data$score <= cutoff, "Low", "High")
    test_data$group <- ifelse(test_data$score <= cutoff, "Low", "High")
    train_data$group <- factor(train_data$group, levels = c("Low", "High"))
    test_data$group <- factor(test_data$group, levels = c("Low", "High"))
    
    tryCatch({
      cox_model <- coxph(Surv(futime, fustat) ~ score + Stage + Grade + New_tumor_event, 
                         data = train_data)
      test_risk <- predict(cox_model, newdata = test_data, type = "risk")
      
      roc <- timeROC(
        T = test_data$futime,
        delta = test_data$fustat,
        marker = test_risk,
        cause = 1,
        times = times,
        iid = TRUE
      )
      
      auc_matrix[valid_count + 1, ] <- roc$AUC
      roc_list[[valid_count + 1]] <- list(FP = roc$FP, TP = roc$TP)
      
      surv_diff <- survdiff(Surv(futime, fustat) ~ group, data = test_data)
      pval <- 1 - pchisq(surv_diff$chisq, df = 1)
      pval_vec[valid_count + 1] <- pval
      
      valid_count <- valid_count + 1
    }, error = function(e) {})
  }
  
  if (valid_count == 0) {
    cat("No valid iterations.\n")
    return(NULL)
  }
  
  auc_matrix <- auc_matrix[1:valid_count, , drop = FALSE]
  pval_vec <- pval_vec[1:valid_count]
  roc_list <- roc_list[1:valid_count]
  
  mean_auc <- colMeans(auc_matrix, na.rm = TRUE)
  sd_auc <- apply(auc_matrix, 2, sd, na.rm = TRUE)
  mean_pval <- mean(pval_vec, na.rm = TRUE)
  
  ## ==== 插值平均ROC曲线 ====
  fixed_fpr <- seq(0, 1, by = 0.01)
  roc_interp_list <- vector("list", length = length(times))
  names(roc_interp_list) <- paste0(times, "-year")
  
  for (j in 1:length(times)) {
    tpr_mat <- matrix(NA, nrow = valid_count, ncol = length(fixed_fpr))
    
    for (i in 1:valid_count) {
      FP_i <- roc_list[[i]]$FP[, j]
      TP_i <- roc_list[[i]]$TP[, j]
      
      if (any(is.na(FP_i))) next
      
      interp_tpr <- approx(FP_i, TP_i, xout = fixed_fpr, ties = mean, rule = 2)$y
      tpr_mat[i, ] <- interp_tpr
    }
    
    roc_interp_list[[j]] <- colMeans(tpr_mat, na.rm = TRUE)
  }
  
  return(list(
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    mean_pval = mean_pval,
    n = valid_count,
    auc_matrix = auc_matrix,
    pval_vec = pval_vec,
    fixed_fpr = fixed_fpr,
    roc_interp_list = roc_interp_list,
    roc_list = roc_list,
    times = times
  ))
}

# ==== 执行Bootstrap验证 ====
result <- bootstrap_validation(data_all, n_boot = 100)

# ==== 绘制带置信区间的ROC曲线 ====
if (!is.null(result)) {
  pdf("Average_ROC_without_CI.pdf", width = 6, height = 6)
  colors <- c("#0066FF", "#FF0000", "#FF9900")
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = "Average ROC Curve")
  abline(0, 1, col = "gray", lty = 2)
  
  # 为每个时间点绘制执行区间
  for (j in 1:length(result$times)) {
    # 提取所有迭代的TPR值
    all_tpr <- matrix(NA, nrow = result$n, ncol = length(result$fixed_fpr))
    
    for (i in 1:result$n) {
      FP_i <- result$roc_list[[i]]$FP[, j]
      TP_i <- result$roc_list[[i]]$TP[, j]
      
      if (any(is.na(FP_i))) next
      
      interp_tpr <- approx(FP_i, TP_i, xout = result$fixed_fpr, ties = mean, rule = 2)$y
      all_tpr[i, ] <- interp_tpr
    }
    
    # 计算95%置信区间
    ci_lower <- apply(all_tpr, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_upper <- apply(all_tpr, 2, quantile, probs = 0.975, na.rm = TRUE)
    

    # 绘制平均ROC曲线
    lines(result$fixed_fpr, result$roc_interp_list[[j]], col = colors[j], lwd = 2)
  }
  
  # 添加图例
  legend_labels <- sprintf("%d-year AUC = %.3f (%.3f-%.3f)", 
                           result$times, 
                           result$mean_auc,
                           result$mean_auc - 1.96*result$sd_auc/sqrt(result$n),
                           result$mean_auc + 1.96*result$sd_auc/sqrt(result$n))
  
  legend("bottomright", legend = legend_labels,
         col = colors, lwd = 2, lty = 1, bty = "n", cex = 0.8,
         title = "Time (95% CI)")
  
  dev.off()
}

# ==== 生存曲线可视化 ====
select_top_survival_plots <- function(data_all, result, cutoff, k = 3) {
  auc_dist <- apply(result$auc_matrix, 1, function(x) sum((x - result$mean_auc)^2))
  selected_ids <- order(auc_dist)[1:k]
  
  for (i in seq_along(selected_ids)) {
    set.seed(selected_ids[i])
    idx <- sample(1:nrow(data_all), size = floor(0.6 * nrow(data_all)))
    train_data <- data_all[idx, ]
    test_data <- data_all[-idx, ]
    
    train_data$group <- ifelse(train_data$score <= cutoff, "Low", "High")
    test_data$group <- ifelse(test_data$score <= cutoff, "Low", "High")
    train_data$group <- factor(train_data$group, levels = c("Low", "High"))
    test_data$group <- factor(test_data$group, levels = c("Low", "High"))
    
    # 训练集生存曲线
    train_fit <- survfit(Surv(futime, fustat) ~ group, data = train_data)
    train_p <- 1 - pchisq(survdiff(Surv(futime, fustat) ~ group, data = train_data)$chisq, df = 1)
    pval_train_text <- ifelse(train_p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", train_p)))
    
    train_plot <- ggsurvplot(train_fit, data = train_data, 
                             pval = pval_train_text, 
                             palette = c("#0066FF", "#FF0000"),
                             legend.title = "Train Group", 
                             risk.table = TRUE, 
                             title = paste0("Training Set - Iter ", selected_ids[i]))
    ggsave(paste0("Train_Surv_Iter", selected_ids[i], ".pdf"), 
           train_plot$plot, width = 7, height = 6)
    
    # 验证集生存曲线
    test_fit <- survfit(Surv(futime, fustat) ~ group, data = test_data)
    test_p <- 1 - pchisq(survdiff(Surv(futime, fustat) ~ group, data = test_data)$chisq, df = 1)
    pval_test_text <- ifelse(test_p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", test_p)))
    
    test_plot <- ggsurvplot(test_fit, data = test_data, 
                            pval = pval_test_text, 
                            palette = c("#0066FF", "#FF0000"),
                            legend.title = "Test Group", 
                            risk.table = TRUE, 
                            title = paste0("Validation Set - Iter ", selected_ids[i]))
    ggsave(paste0("Test_Surv_Iter", selected_ids[i], ".pdf"), 
           test_plot$plot, width = 7, height = 6)
  }
}

# 执行生存曲线可视化
if (!is.null(result)) {
  select_top_survival_plots(data_all, result, cutoff = cutoff, k = 3)
}

# ==== 结果汇总输出 ====
if (!is.null(result)) {
  cat("\n=== Final Results ===\n")
  cat("Valid bootstrap iterations:", result$n, "\n")
  for (j in 1:length(result$times)) {
    cat(sprintf("%d-year AUC: %.3f (95%% CI: %.3f-%.3f)\n", 
                result$times[j], 
                result$mean_auc[j],
                result$mean_auc[j] - 1.96*result$sd_auc[j]/sqrt(result$n),
                result$mean_auc[j] + 1.96*result$sd_auc[j]/sqrt(result$n)))
  }
  cat("Mean log-rank p-value:", round(result$mean_pval, 4), "\n")
}