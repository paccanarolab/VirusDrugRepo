library(Polychrome)

source("Data collection/load_additional_biological_data.R")
source("Data analysis/biological_interpretation.R")

D <- NULL
loadEmbedding()
loadBiologicalData(sel_drugs = rownames(D), sel_viruses=rownames(V), 
                   sel_proteins = colnames(P), sel_genes=colnames(U))
computeSimilarities()

# Plot showing correlation between drug embedding similarity and chemical
# similarity (Fig 4e)
pdf(paste0("Data analysis/Plots/", "cor_chem.pdf"))
p <- plotCorSim(simD, chem_sim, label="Chemical similarity")
dev.off()
print(p)

# Plot showing correlation between protein embedding similarity and semantic
# similarity (GO Biological Processes) (Fig 4f)
pdf(paste0("Data analysis/Plots/", "cor_bp.pdf"))
p <- plotCorSim(simP[targets, targets], entrez_bp, label="Semantic similarity")
dev.off()
print(p)

# Plot showing correlation between protein embedding similarity and semantic
# similarity (GO Biological Processes) (Supplementaty Fig S8a)
pdf(paste0("Data analysis/Plots/", "cor_cc.pdf"))
p <- plotCorSim(simP[targets, targets], entrez_cc, label="Semantic similarity")
dev.off()
print(p)

# Plot showing correlation between protein embedding similarity and semantic
# similarity (GO Biological Processes) (Supplementaty Fig S8b)
pdf(paste0("Data analysis/Plots/", "cor_mf.pdf"))
p <- plotCorSim(simP[targets, targets], entrez_mf, label="Semantic similarity")
dev.off()
print(p)

# Heatmap of average drug embedding similarities inter and intra ATC categories
pdf(paste0("Data analysis/Plots/", "heatmap_drug_sim_ATC.pdf"), width=8)
p <- plotSim(avgSimD, low="white", mid=NULL)
dev.off()
print(paste("pvalue:", p))

# Barplots comparing embedding similarities for drug-virus pairs, protein pairs
# drug pairs, and gene pairs between pairs that do not share pathways and
# pairs that share pathways (Fig 4a)
pdf(paste0("Data analysis/Plots/", "barplots_similarities_pathways.pdf"), width=10)
barPlots("No shared\npathways\n", "Shared\npathways")
dev.off()

# Virus plots ------------------------------------------------------------------

# Compute cosine distance between virus embeddings 
dist1 <- dist(1-cos_similarities_norm_V[virus_3, virus_3])
hmethod <- "single"
h1 <- hclust(dist1, method=hmethod)
# Create virus dendogram
d1 <- as.dendrogram(h1)

colors <- rainbow(ncol(virus_family))
colors <- glasbey.colors(ncol(virus_family)+1)
colors[1] <- colors[length(colors)]
colors <- colors[-length(colors)]

names(colors) <- colnames(virus_family)
colors["Picornaviridae"] <- "#dc7a00"

colLab <- function(n) {
  if(is.leaf(n))
  {
    a <- attributes(n)
    family <- group_per_virus[[a$label]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colors[family], lab.cex=.7, 
                                            col=colors[family], cex=.7, pch=16 ))
  }
  return(n)
}

# Select families with at least 3 viruses
sel <- which(colSums(virus_family[virus_3,]) > 0)
names(colors) <- colnames(virus_family[virus_3,sel])

# Finally I just have to apply this to my dendrogram
dL <- dendrapply(d1, colLab)

# Plot dendogram of viruses (Fig 4d)
pdf(paste0("Data analysis/Plots/", "dendogram.pdf"), width=15)
plot(dL)
dev.off()

# Plot 3d representation of virus signatures (Fig 4c)
labels <- unlist(group_per_virus)
labels <- labels[rownames(V)]
labels["SFV"] <- "Retroviridae"
pdf(paste0("Data analysis/Plots/", "viruses_tsne_3dim_fam.pdf"), width=8)
plotComponents(V[virus_3, ], labels[virus_3], palete=colors[unique(labels[virus_3])])
dev.off()



