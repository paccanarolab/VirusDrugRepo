
source("Data collection/load_evaluation_data.R")
source("Data analysis/evaluation.R")

# Loading prediction scores ----------------------------------------------------

# loading prediction scores with different lambda values
loadPredctionsLambdas()

# Computing AUC and recall -----------------------------------------------------

# Select viruses with evaluation data available
loadDrugVirusAssociations()
i <- which(colSums(drug_virus) > 0)
sel_viruses <- colnames(drug_virus)[i]
drugs <- rownames(our_method)


# Compute AUC and recall (Fig S9)
eval_lambdas = evaluateMethods(scores_lambdas, sel_viruses, 
                                   rownames(scores_lambdas_diff[[1]]))

# Compute AUC and recall (Fig S10)
eval_lambdas_diff = evaluateMethods(scores_lambdas_diff, sel_viruses, 
                                        rownames(scores_lambdas_diff[[1]]))

# Generate plots ---------------------------------------------------------------

colors <- c("#77a451", "#9765ca","#c7733b", "#6c93c5","#c65b80")
results_folder <- "Data analysis/Plots/"

# Supplementary Fig S9
pdf(paste0(results_folder, "AUC_lambdas.pdf"), width = 8)
boxplot(eval_lambdas$auc, pch=20, col=colors, ylab="AUC", cex.lab=1.5, 
        cex.axis=1.5)
dev.off()

# Supplementary Fig S10
pdf(paste0(results_folder, "AUC_lambdas_diff.pdf"), width = 8)
boxplot(eval_lambdas_diff$auc, pch=20, col=colors, ylab="AUC", cex.lab=1.5, 
        cex.axis=1.5)
dev.off()

