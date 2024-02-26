source("Data analysis/cosine_similarity.R")

# Load drug-target associations
load("Data collection/Data/drug_targets_Gysi.RData")
# Load virus-host protein associations
load("Data collection/Data/virus_host_Gysi.RData")

# Load embedding learned by the model
loadEmbedding <- function() {

  load("Data analysis/Embedding/selectedEmbedding.RData")
  
  V <<- results$V
  
  D <<- results$D
  
  P <<- results$P
  
  U <<- results$U
}

# Compute average cosine similarity per group
computeAvgSimilarities <- function(groups, data, fun=cosSimilarity) {
  similarities <- fun(data)
  diag(similarities) <- NA
  avgSim <- matrix(NA, length(groups), length(groups))
  rownames(avgSim) <- names(groups)
  colnames(avgSim) <- names(groups)
  for (i1 in 1:(length(groups))) {
    g1 <- names(groups)[i1]
    sel1 <- groups[[g1]]
    for (i2 in i1:length(groups)){
      g2 <- names(groups)[i2]
      sel2 <- groups[[g2]]
      s <- mean(similarities[sel1, sel2], na.rm=TRUE)
      avgSim[g1, g2] <- s
    }
  }
  return(avgSim)
}

# Compute drug, virus, protein, and gene similarities
computeSimilarities <- function() {
  simD <<- cosSimilarity(D)
  simV <<- cosSimilarity(V)
  simP <<- cosSimilarity(t(P))
  simU <<- cosSimilarity(t(U))
  simDV <<- cosSimilarityAB(D, V)
  avgSimD <<- computeAvgSimilarities(drug_groups, D)
  avgSimV <<- computeAvgSimilarities(virus_groups, V)
}

# Build an association matrix indicating whether two items share at least one
# group
buildAssociationMatrix <- function(groups, filter=NULL) {
  items <- unique(unlist(groups))
  A <- matrix(0, length(items), length(groups))
  rownames(A) <- items
  colnames(A) <- names(groups)
  for (g in names(groups)) {
    sel <- groups[[g]]
    if (!is.null(filter))
      sel <- intersect(sel, filter)
    if (length(sel) > 0)
      A[sel, g] <- 1
  }
  return(A)
}

# If converts a list of groups of elements into a list of elements indicating 
# which groups they belong to.
groupPerElement <- function(groups) {
  group_per_element <- list()
  for (g in names(groups)) {
    sel_elements <- groups[[g]]
    for (element in sel_elements)
      group_per_element[[element]] <- c(group_per_element[[element]], g)
  }
  return(group_per_element)
}

# Load data based on aditional biological information
loadBiologicalData <- function(sel_drugs=NULL, sel_viruses=NULL, 
                               sel_proteins=NULL, sel_genes=NULL) {
  
  dir <- "Data collection/Data/"
  file <- "DrugBank_ATC.tsv" # ATC categories obtained from drugbank
  
  # Drugs ---------------------------------------------------------------------
  
  atc_data <- as.matrix(read.csv(paste0(dir, file), sep="\t"))
  
  if (!is.null(sel_drugs)) {
    i <- which(atc_data[,1] %in% sel_drugs)
    atc_data <- atc_data[i, ]
  }
  
  atcs <- unique(atc_data[,2])
  
  text <- c()
  
  sizes <- c()
  
  getFirstLetter <- function(word, end=1) {
    return(substr(word, 1, 1))
  }
  
  firstLetters <- unique(unlist(lapply(atcs, getFirstLetter)))
  drug_groups <- list()
  for (atc in firstLetters) {
    i <- which(startsWith(atc_data[,2], atc))
    drugs <- atc_data[i,1]
    drug_groups[[atc]] <- unique(drugs)
  }
  drug_groups <<- drug_groups
  drug_atc <<- buildAssociationMatrix(drug_groups, filter=sel_drugs)
  
  load("Data collection/Data/drug_kegg_pathways.RData")
  drug_kegg <<- buildAssociationMatrix(drug_pathways_adj)
  
  load("Data collection/Data/drug_chem_sim.RData")
  chem_sim <<-chem_sim 
  
  load("Data collection/Data/drug_pathways.RData")
  drug_drugbank <<- drug_pathways_matrix
  drug_drugbank_sim <<- drug_pathways_matrix%*%t(drug_pathways_matrix)
  
  
  # Virus pathways ---------------------------------------------------------------
  
  # KEGG pathways that are differentially expressed in viral infections. 
  # We selected pathways that are differentially expressed for at leat 3 
  # different viruses.
  load("Data collection/Data/virus_diff_expressed_pathways.RData")
  
  sizes_virus_pathways <- unlist(lapply(virus_pathways_adj, length))
  virus_pathways_adj <- virus_pathways_adj[which(sizes_virus_pathways >= 3)]

  virus_kegg <<- buildAssociationMatrix(virus_pathways_adj, filter=sel_viruses)
  virus_kegg_sim <<- virus_kegg%*%t(virus_kegg)
  
  # Virus family ---------------------------------------------------------------
  load("Data collection/Data/virus_families.RData")
  virus_groups <<- virus_groups # Virus families
  sizes_virus_groups <- unlist(lapply(virus_groups, length))
  i <- which(sizes_virus_groups >= 3)
  virus_groups_3 <<- virus_groups[i] # Virus families with at leat 3 viruses
  virus_3 <<- intersect(rownames(V), unlist(virus_groups_3))
  virus_family <<- buildAssociationMatrix(virus_groups, filter=sel_viruses)
  group_per_virus <<- groupPerElement(virus_groups)
  
  # Virus embedding similarities
  cos_similarities_norm_V <<- cosSimilarity(V, norm = TRUE)
  
  # Protein pathways -------------------------------------------------------------
  load("Data collection/Data/kegg_pathways.RData")
  targets <<- intersect(colnames(drug_targets), colnames(P))
  host <<- intersect(colnames(virus_host), colnames(P))
  
  protein_kegg_targets <<- buildAssociationMatrix(kegg_pathways_entrez, 
                                                  filter=targets)
  protein_kegg_targets <<- protein_kegg_targets[intersect(targets, 
                                              rownames(protein_kegg_targets)), ]
  protein_kegg_targets_sim <<- protein_kegg_targets%*%t(protein_kegg_targets)
  protein_kegg <<- buildAssociationMatrix(kegg_pathways_entrez)
  protein_kegg_sim <<- protein_kegg%*%t(protein_kegg)
  gene_kegg <<- buildAssociationMatrix(kegg_pathways_symbol)
  gene_kegg_sim <<- gene_kegg%*%t(gene_kegg)
  
  # Protein semantic similarities ----------------------------------------------
  # Load semantic similarities between proteins according to GO biolocial 
  # processes
  load("Data collection/Data/sim_bp_entrez.RData")
  # Load semantic similarities between proteins according to GO cellular  
  # components
  load("Data collection/Data/sim_cc_entrez.RData")
  # Load semantic similarities between proteins according to GO molecular 
  # function
  load("Data collection/Data/sim_mf_entrez.RData")
  
  entrez_bp <<- entrez_bp
  entrez_cc <<- entrez_cc
  entrez_mf <<- entrez_mf
  
  # drug virus -----------------------------------------------------------------
  load("Data collection/Data/drug_kegg_pathways.RData")
  drug_kegg_adj <<- buildAssociationMatrix(drug_pathways_adj, 
                                           filter=rownames(D))
  
  load("Data collection/Data/pvalues_drug_virus_kegg.RData")
  p_adj <- p.adjust(pvalues_drug_virus_kegg, method="fdr")
  sel <-  which(p_adj < 0.05)
  simDV_kegg <<- drug_kegg[,names(sel)]%*%t(virus_kegg[,names(sel)])
}
