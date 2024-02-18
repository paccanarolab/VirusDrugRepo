
source("Data collection/load_evaluation_data.R")
source("Data analysis/evaluation.R")

# Computing AUC and recall -----------------------------------------------------

dir_expr <- "Data analysis/Embedding/Different initializations - alpha = 1/"
dir_noexpr <- "Data analysis/Embedding/Different initializations - alpha = 0/"
dir_assess <- "Data analysis/Embedding/Different initializations - Assessment set/"

files_expr <- dir(dir_expr)
files_noexpr <- dir(dir_noexpr)
files_assess <- dir(dir_assess)

# Compute AUC and Recall for experiments without gene expression data (alpha = 0)
eval_noexpr <- evaluateFiles(dir_noexpr, files_noexpr)

# Compute AUC and Recall for experiments with gene expression data (alpha = 1)
eval_expr <- evaluateFiles(dir_expr, files_expr)

# Compute AUC and Recall for experiments with only viruses from the assessment
# set
eval_assess <- evaluateFiles(dir_assess, files_assess)

# Generate plots ---------------------------------------------------------------

colors <- c("#77a451", "#9765ca","#c7733b", "#6c93c5","#c65b80")
labels <- c("Without gene expression",  "With gene expression")
results_folder <- "Data analysis/Plots/"

boxplotCompare <- function(x, y, ylab="", labels=NULL) {
  data <- list(x, y)
  names(data) <- labels
  boxplot(data, col=colors[c(4,3)], ylab=ylab, cex.axis=1.5, cex.lab=1.5)
}

## Plots Figure 3 - effect of gene expression
# Figure 3a
pdf(paste0(results_folder, "Expr_median_AUC.pdf"), width=8)
par(mar=c(3,5,3,1)+.1)
boxplotCompare(apply(eval_noexpr$auc, 2, median), apply(eval_expr$auc, 2, 
                                                        median), 
               
               ylab="Median AUC", labels=labels)
dev.off()

# Figure 3b
pdf(paste0(results_folder, "Expr_median_Recall.pdf"), width=8)
par(mar=c(3,5,3,1)+.1)
boxplotCompare(apply(eval_noexpr$recall, 2, median), apply(eval_expr$recall, 2, 
                                                           median), 
               ylab="Median recall", labels=labels)
dev.off()

## Plots Supplementary Fig S3 - effect of adding viruses

labels <- c("Assessment set",  "Original set")

# Fig S3a
pdf(paste0(results_folder, "assess_median_AUC.pdf"), width=8)
par(mar=c(3,5,3,1)+.1)
boxplotCompare(apply(eval_assess$auc, 2, median), apply(eval_expr$auc, 2, 
                                                          median), 
               ylab="Median AUC", labels=labels)
dev.off()

# Fig S3b
pdf(paste0(results_folder, "assess_median_Recall.pdf"), width=8)
par(mar=c(3,5,3,1)+.1)
boxplotCompare(apply(eval_assess$recall, 2, median), apply(eval_expr$recall, 
                                                             2, median), 
               ylab="Median recall", labels=labels)
dev.off()

