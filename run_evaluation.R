
source("Data collection/load_evaluation_data.R")
source("Data analysis/evaluation.R")

# Loading prediction scores ----------------------------------------------------

our_method <- getPredictions("Data analysis/Embedding/selectedEmbedding.RData")

loadPredictionsCompetitors()

scores_main <- list("Standard network medicine approach"=santos,
                 "Our method"=our_method)
scores_sup_li <- list("Li et al"=li_et_al, "Our method"=our_method)

scores_sup_saverunner <- list("SaveRUNNER"=save_runner, "Our method"=our_method)


# Computing AUC and recall -----------------------------------------------------

# Select viruses with evaluation data available
loadDrugVirusAssociations()
i <- which(colSums(drug_virus) > 0)
sel_viruses <- colnames(drug_virus)[i]
drugs <- rownames(our_method)

# Compute AUC and recall (Fig 2 and Fig S)
eval_main <- evaluateMethods(scores_main, sel_viruses, drugs)

# Select viruses with predictions by Li et al.
i <- which(colSums(li_et_al, na.rm = TRUE) > 0)
sel_viruses_li <- colnames(li_et_al)[i]

# Select drugs with predictions by Li et al.
i <- which(!is.na(li_et_al[, "HCV"]))
sel_drugs_li <- rownames(li_et_al)[i]

# Compute AUC and recall (Fig S2)
eval_sup_li <- evaluateMethods(scores_sup_li, sel_viruses_li, sel_drugs_li)

# Compute AUC and recall (Fig SX)
eval_sup_saverunner <- evaluateMethods(scores_sup_saverunner, sel_viruses, drugs)

# Generate plots ---------------------------------------------------------------

colors <- c("#77a451", "#9765ca","#c7733b", "#6c93c5","#c65b80")

results_folder <- "Data analysis/Plots/"

# Fig 2a
pdf(paste0(results_folder, "AUC_main.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$auc, pch=20, col=colors[2:3], ylab="AUC", cex.lab=1.5, 
        cex.axis=1.5, names=c("Standard network\nmedicine approach","Our method"))
dev.off()


# Fig 2b
pdf(paste0(results_folder, "Recall_main.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$recall, pch=20, col=colors[2:3], ylab="Recall@150", cex.lab=1.5, 
        cex.axis=1.5, names=c("Standard network\nmedicine approach","Our method"))
dev.off()

# Supplementary Fig S1a
pdf(paste0(results_folder, "Recall_20.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$recall_20, pch=20, col=colors[2:3], ylab="Recall@20", cex.lab=1.5, 
        cex.axis=1.5, names=c("Standard network\nmedicine approach","Our method"))
dev.off()

# Supplementary Fig S1b
pdf(paste0(results_folder, "Recall_50.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$recall_50, pch=20, col=colors[2:3], ylab="Recall@50", cex.lab=1.5, 
        cex.axis=1.5, names=c("Standard network\nmedicine approach","Our method"))
dev.off()

# Supplementary Fig S1c
pdf(paste0(results_folder, "Recall_100.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$recall_100, pch=20, col=colors[2:3], ylab="Recall@100", 
        cex.lab=1.5, cex.axis=1.5, names=c("Standard network\nmedicine approach",
                                           "Our method"))
dev.off()

# Supplementary Fig S1d
pdf(paste0(results_folder, "Recall_200.pdf"))
par(mar=c(3,5,3,1)+.1, mgp=c(3,2,0))
boxplot(eval_main$recall_200, pch=20, col=colors[2:3], ylab="Recall@200", 
        cex.lab=1.5, cex.axis=1.5, names=c("Standard network\nmedicine approach",
                                           "Our method"))
dev.off()

# Supplementary Fig S4a
pdf(paste0(results_folder, "AUC_sup_Li.pdf"), width=8)
boxplot(eval_sup_li$auc, pch=20, col=c(colors[5], colors[3]), ylab="AUC", 
        cex.lab=1.5, cex.axis=1.5, names=c("Li at al", "Our method"))
dev.off()

# Supplementary Fig S4b
pdf(paste0(results_folder, "Recall_sup_Li.pdf"), width=8)
boxplot(eval_sup_li$recall, pch=20, col=c(colors[5], colors[3]), 
        ylab="Recall@150", cex.lab=1.5, cex.axis=1.5, names=c("Li at al", 
                                                              "Our method"))
dev.off()


# Supplementary Fig S3a
pdf(paste0(results_folder, "AUC_sup_SaveRUNNER.pdf"), width=8)
boxplot(eval_sup_saverunner$auc, pch=20, col=c(colors[1], colors[3]), ylab="AUC", 
        cex.lab=1.5, cex.axis=1.5, names=c("SaveRUNNER", "Our method"))
dev.off()

# Supplementary Fig S3b
pdf(paste0(results_folder, "Recall_sup_SaveRUNNER.pdf"), width=8)
boxplot(eval_sup_saverunner$recall, pch=20, col=c(colors[1], colors[3]), 
        ylab="Recall@150", cex.lab=1.5, cex.axis=1.5, names=c("SaveRUNNER", 
                                                              "Our method"))
dev.off()
