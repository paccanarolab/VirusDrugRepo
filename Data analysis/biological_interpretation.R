require(ggplot2)
require(reshape2)
require(Rtsne)
require(scatterplot3d)

# Compare embedding similarity between items that have similarity bigger than
# zero in "sim"
compareSim <- function(embedding_sim, sim, 
                       value=0, filename="", symmetric=TRUE) {
  
  inter1 <- intersect(rownames(embedding_sim), rownames(sim))
  inter2 <- intersect(colnames(embedding_sim), colnames(sim))
  
  embedding_sim <- embedding_sim[inter1, inter2]
  sim <- sim[inter1, inter2]
  
  if (symmetric) {
    embedding_sim <- embedding_sim[upper.tri(embedding_sim)]
    sim <- sim[upper.tri(sim)]
  }
  
  sel <- which(sim > value)
  
  pval <- wilcox.test(embedding_sim[-sel], embedding_sim[sel], alternative="less")$p.value
  {if (pval < 0.0005)
    res <- " < 0.0005"
    else
      res <- paste0(" = ",  round(pval, 3))
  }
  
  return(list("pvalue"=pval, "sim1"=embedding_sim[-sel], "sim2"=embedding_sim[sel]))
}

# Compute the mean and confidence interval around the mean
meanCi <- function(comp, legend1="0", legend2=">0") {
  means <- unlist(lapply(comp, mean))[-1]
  names(means) <- c(legend1, legend2)
  sds <- unlist(lapply(comp, sd))[-1]
  ns <- unlist(lapply(comp, length))[-1]
  ci <- sds*1.96/sqrt(ns)
  return(list("mean"=means, "ci"=ci))
}

# barPlots used for building Fig 4a
barPlots <- function(legend1="0", legend2=">0") {
  comp_drug_drugbank <- compareSim(simD, drug_drugbank_sim)
  comp_prot_kegg <- compareSim(simP, protein_kegg_sim)
  comp_gene_kegg <- compareSim(simU, gene_kegg_sim)
  comp_drug_virus_kegg <- compareSim(simDV, simDV_kegg, symmetric = FALSE)
  stat_drug_virus_kegg <- meanCi(comp_drug_virus_kegg)
  stat_drug_drugbank <- meanCi(comp_drug_drugbank)
  stat_protein_kegg <- meanCi(comp_prot_kegg)
  stat_gene_kegg <- meanCi(comp_gene_kegg)
  means <- matrix(NA, 2, 4)
  rownames(means) <- c(legend1, legend2)
  colnames(means) <- c("Drug-virus pairs", "Drug pairs", "Protein pairs", 
                       "Gene pairs")
  ci <- means
  means[, "Drug-virus pairs"] <- stat_drug_virus_kegg$mean
  ci[, "Drug-virus pairs"] <- stat_drug_virus_kegg$ci
  means[, "Drug pairs"] <- stat_drug_drugbank$mean
  ci[, "Drug pairs"] <- stat_drug_drugbank$ci
  means[, "Protein pairs"] <- stat_protein_kegg$mean
  ci[, "Protein pairs"] <- stat_protein_kegg$ci
  means[, "Gene pairs"] <- stat_gene_kegg$mean
  ci[, "Gene pairs"] <- stat_gene_kegg$ci
  
  lim <- max(means + ci + 0.1)
  palete1 <- c("#a670b0","#a58d48")
  palete2 <- c("#f79891","#58deff")
  palete3 <- c("#6c93c5", "#c7733b")
  par(mar=c(5,6,4,1)+.1)
  fig <- barplot(means , beside=T , legend.text=T,col=palete3, ylim=c(0,lim) , 
                 ylab="Similarity", cex.axis = 1.5, cex.names = 1.5, 
                 cex.lab=1.5, main="Average cosine similarities between embeddings")
  arrows(fig,means-ci, fig, means+ci, angle=90, code=3, length=0.1, lwd = 2)
}

# Scatter plot and correlation test between similarities. It was used for
# Fig 4e,f and Fig S8a,b
plotCorSim <- function(embedding_sim, sim, label="", n=100) {
  inter <- intersect(rownames(embedding_sim), rownames(sim))
  embedding_sim <- embedding_sim[inter, inter]
  sim <- sim[inter, inter]
  
  embedding_sim <- embedding_sim[upper.tri(embedding_sim)]
  sim <- sim[upper.tri(sim)]
  
  ntotal <- length(sim)
  size <- floor(ntotal/n)
  
  o <- order(sim)
  x <- c()
  y <- c()
  for (i in 1:n) {
    j <- (i-1)*size + 1
    x[i] <- mean(sim[o[j:(j+size-1)]])
    y[i] <- mean(embedding_sim[o[j:(j+size-1)]])
  }
  
  plot(x, y, xlab=label, ylab="Embedding similarity", cex.axis=1.5, cex=2, cex.lab=1.5, pch=19)
  fit <- lm(y ~ x)
  lines(x, fit$fitted.values, col="#6c93c5", lwd=5)
  return(cor.test(x, y))
}

# Build heatmap for Fig 4b
plotSim <- function(avgSim, show_pvalue=FALSE, norm_col=TRUE, show_mean=FALSE,
                    low="blue", mid="white", high="red") {
  if (norm_col) {
    diags <- diag(avgSim)
    avgSim[which(is.na(avgSim))] <- 0
    avgSim <- avgSim + t(avgSim)
    diag(avgSim) <- diags
    mins <- apply(avgSim, 2, function(x){return(min(x, na.rm=TRUE))})
    # mins <- 0
    avgNorm <- t(avgSim) - mins
    maxs <- apply(avgNorm, 1, function(x){return(max(x, na.rm=TRUE))})
    avgSim <-  t(avgNorm/maxs) 
    #avgSim[lower.tri(avgSim)] <- NA
  }
  
  o <- order(diag(avgSim))
  avgSim <- avgSim[o,o]
  
  melted_cormat <- melt(avgSim[,ncol(avgSim):1], na.rm = TRUE)
  # Heatmap
  minv <- min(avgSim, na.rm = TRUE)
  maxv <- max(avgSim, na.rm=TRUE)
  meanv <- max(avgSim, na.rm=TRUE)
  avg_remaining <- avgSim[upper.tri(avgSim)]
  avg_diag <- diag(avgSim)
  
  pvalue <- wilcox.test(avg_diag, avg_remaining, alternative="greater")$p.value
  tmp <- ""
  if (show_pvalue)
    tmp <- paste0("\nP-value: ", round(pvalue, 4))
  title <- paste0("Mean similarity: ", round(mean(avg_remaining), 3), " (lower tri), ", 
                  round(mean(avg_diag), 3), " (diag)", tmp)
  if (!show_mean) {
    if (norm_col)
      title <- "Average similarities normalized by column"
    else
      title <- "Average similarities"
  }
  g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low =low, high =high, mid=mid, midpoint = minv + (maxv-minv)/2,
                         limit = c(minv,maxv), space = "Lab", 
                         name="Cosine similarity") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 16, hjust = 1), 
          axis.text.y=element_text(size=16), plot.title=element_text(size=16), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16), axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    coord_fixed() +
    ggtitle(title)
  print(g)
  return(pvalue)
}

# For each KEGG pathway that is differentially expressed in cell lines infected 
# by viruses and cell lines treated with drugs, it tests whether viruses and 
# drugs are closer in the embedding when they share the pathways.
# Then, the p-values are adjusted for multiple testing by FDR.
# We found 3 pathways with significant adjusted p-value.
selectDrugVirusKeggPathways <- function() {
  pvalues_drug_virus_kegg <- c()
  for (pathway in pathways) {
    sim <- drug_kegg_adj[,pathway]%*%t(virus_kegg[, pathway])
    rownames(sim) <- rownames(drug_kegg_adj)
    if (sum(sim) > 1) {
      p <- compareSim(simDV, sim, symmetric = FALSE)$pvalue
      pvalues_drug_virus_kegg[pathway] <- p
      print(paste(pathway, p))
    }
  }
  p_adj <- p.adjust(pvalues_drug_virus_kegg, method="fdr")
  sel <-  which(p_adj < 0.05)
  save(pvalues_DV_kegg_adj, file=paste0("Data collection/Data/pvalues_drug_virus_kegg.RData"))
  return(names(sel))
}

# Dimension reduction by tsne (Fig 4c)
plotComponents <- function(x, labels, palete=NULL, fam_text=NULL) {
  set.seed(0)
  dim = 3
  tsne <- Rtsne(x, dim=dim)
  x1 <- tsne$Y[, 1]
  x2 <- tsne$Y[, 2]
  x3 <- tsne$Y[, 3]
  data <- data.frame("comp1"=x1, "comp2"=x2, "comp3" = x3)
  data <- cbind(data, "Label"=labels)
  colors <- palete
  names(colors) <- unique(labels)
  cols <- colors[labels]
  s3d <- scatterplot3d(data[,1:3],  color=cols, pch=20, cex.symbols = 2.3)
  legend("topright", legend = unique(labels),
         col =  unique(cols), pch=20, pt.cex=2.3, bg="white")
  if (!is.null(fam_text)) {
    scatter3D(x1, x2, x3, col = cols, add = FALSE, pt.cex=2.3)
    i <- which(rownames(x) %in% virus_groups[[fam_text]])
    text3D(x1[i], 
           x2[i],  x3[i],           
           labels = rownames(x)[i],               
           cex = 1, 
           pos = 4)  
  }
  
}

