
source("Data collection/load_evaluation_data.R")
source("Data analysis/evaluation.R")

# Loading prediction scores ----------------------------------------------------

our_method <- getPredictions("Data analysis/Embedding/selectedEmbedding.RData")

loadPredictionsCompetitors()

scores_main <- list("SaveRUNNER"=save_runner, "Santos"=santos,
                 "Our method"=our_method)
scores_sup <- list("Li et al"=li_et_al, "Our method"=our_method)


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
eval_sup <- evaluateMethods(scores_sup, sel_viruses_li, sel_drugs_li)


# Generate plots ---------------------------------------------------------------

colors <- c("#77a451", "#9765ca","#c7733b", "#6c93c5","#c65b80")

results_folder <- "Data analysis/Plots/"

# Fig 2a
pdf(paste0(results_folder, "AUC_main.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$auc, pch=20, col=colors, ylab="AUC", cex.lab=1.5, 
        cex.axis=1.5, names=c("SAveRUNNER","Santos et al","Our method"))
dev.off()

# Fig 2b
pdf(paste0(results_folder, "Recall_main.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$recall, pch=20, col=colors, ylab="Recall@150", cex.lab=1.5, 
        cex.axis=1.5, names=c("SAveRUNNER","Santos et al","Our method"))
dev.off()

# Supplementary Fig S1a
pdf(paste0(results_folder, "Recall_20.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$recall_20, pch=20, col=colors, ylab="Recall@20", cex.lab=1.5, 
        cex.axis=1.5, names=c("SAveRUNNER","Santos et al","Our method"))
dev.off()

# Supplementary Fig S1b
pdf(paste0(results_folder, "Recall_50.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$recall_50, pch=20, col=colors, ylab="Recall@50", cex.lab=1.5, 
        cex.axis=1.5, names=c("SAveRUNNER","Santos et al","Our method"))
dev.off()

# Supplementary Fig S1c
pdf(paste0(results_folder, "Recall_100.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$recall_100, pch=20, col=colors, ylab="Recall@100", 
        cex.lab=1.5, cex.axis=1.5, names=c("SAveRUNNER","Santos et al",
                                           "Our method"))
dev.off()

# Supplementary Fig S1d
pdf(paste0(results_folder, "Recall_200.pdf"))
par(mar=c(3,5,3,1)+.1)
boxplot(eval_main$recall_200, pch=20, col=colors, ylab="Recall@200", 
        cex.lab=1.5, cex.axis=1.5, names=c("SAveRUNNER","Santos et al",
                                           "Our method"))
dev.off()

# Supplementary Fig S2a
pdf(paste0(results_folder, "AUC_sup.pdf"), width=8)
boxplot(eval_sup$auc, pch=20, col=c(colors[5], colors[3]), ylab="AUC", 
        cex.lab=1.5, cex.axis=1.5, names=c("Li at al", "Our method"))
dev.off()

# Supplementary Fig S2b
pdf(paste0(results_folder, "Recall_sup.pdf"), width=8)
boxplot(eval_sup$recall, pch=20, col=c(colors[5], colors[3]), 
        ylab="Recall@150", cex.lab=1.5, cex.axis=1.5, names=c("Li at al", 
                                                              "Our method"))
dev.off()


