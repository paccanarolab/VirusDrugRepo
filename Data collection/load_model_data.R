require(org.Hs.eg.db) # library for obtaining valid gene symbols

# Loading kernel matrix K
# For obtaining matrix K we used the PPI from 
# https://doi.org/10.1073/pnas.2025581118 to build an adjacency matrix. Then
# we used the pStepKernel function from diffStats R package with parameter 
# p = 2 to compute the kenel.
load("~/Data/Gradient Descent/kernel.RData")

# Loading drug target associations obtained originally from 
# https://doi.org/10.1073/pnas.2025581118  
load("~/Dropbox/Projetos/Métodos/Drug repositioning/Paper/code/Data/drug_targets_Gysi.RData")

# Loading virus-host associations obtained from HVDB. Only proteins available
# from PPI at https://doi.org/10.1073/pnas.2025581118 were included
load("~/Dropbox/Projetos/Métodos/Drug repositioning/Paper/code/Data/virus_host_Gysi.RData")

# Getting protein entrez IDs
proteins <- rownames(K)

# Getting Drugbank drug IDs
drugs <- rownames(drug_targets)

# Preprocessing for building matrix A ------------------------------------------

drug_host <- matrix(0, length(drugs), length(proteins))
rownames(drug_host) <- drugs
colnames(drug_host) <- proteins

all_targets <- colnames(drug_targets)

for (drug in drugs) {
  drug_host[drug,] <- apply(drug_targets[drug,all_targets]*K[all_targets, 
                                                             proteins], 2, max)
}

# Preprocessing for building matrix B ------------------------------------------

virus_target <- matrix(0, nrow(virus_host), length(proteins))
rownames(virus_target) <- rownames(virus_host)
colnames(virus_target) <- proteins

host_prots <- colnames(virus_host)

viruses <- rownames(virus_host)
for (virus in viruses) {
  virus_target[virus,] <- apply(virus_host[virus,host_prots]*K[host_prots, 
                                                               proteins],2, max)
}

# Making kernel scores binary for building A and B -----------------------------

pth <- 0.9
th <- quantile(drug_host, pth)
i <- which(drug_host >= th) # selecting 10% highest values
maximum <- max(drug_host, virus_target)
drug_host[i] <- 1
j <- setdiff(1:length(drug_host),i)
drug_host[j] <- 0
th <- quantile(virus_target, pth)
i <- which(virus_target >= th)
virus_target[i] <- 1
j <- setdiff(1:length(virus_target),i)
virus_target[j] <- 0
VH <- virus_target
DH <- drug_host
i <- which(colSums(VH) == 0)
j <- which(colSums(DH) == 0)
sel_proteins <- setdiff(colnames(VH), union(names(i), names(j)))

# Building matrices A and B
A <- DH[, sel_proteins]
B <- VH[, sel_proteins]

# Building matrix G ------------------------------------------------------------

# This function builds a matrix containing z-scores that measure the 
# difference of expression between cell lines infected by a virus and 
# non-infected cell lines. 
# INPUT: a list of viruses
# OUTPUT: a matrix G, with viruses on the rows and genes on the columns. The 
# position G[i,j] indicates the z-score obtained for virus i and gene j, that
# is, the differential expression of gene j during an infection by virus i.
geneExpressionMatrix <- function(viruses) {
  expr_dir <- "~/Dropbox/Projetos/Dados/Viruses expression/Differential expression/"
  G <- NULL
  for (virus in viruses) {  
    
    expr_file <- paste0("zscore_", virus, ".RData")
    file <- paste0(expr_dir, expr_file)
    if (file.exists(file)) {
      load(file)
      
      # Filter out invalid z-scores (nans, NAs, and infite)
      weights <- weights[which(!is.nan(weights))]
      weights <- weights[which(!is.na(weights))]
      weights <- weights[which(!is.infinite(weights))]
      genes <- unique(names(weights))
      genes <- genes[genes != ""]
      expr <- NULL
      for (g in genes) {
        # If more than one probe set was mapped to a sane gene symbol we select
        # the highest z-score
        i <- which(as.character(names(weights)) %in% g)
        if (length(i) > 0) {
          expr[g] <- max(weights[i], na.rm=TRUE)
        }
      }
      # Build matrix G when the first virus is selected
      if (is.null(G)) {
        
        G <- matrix(expr, nrow=1, ncol=length(expr))
        colnames(G) <- names(expr)
        rownames(G) <- virus
      }
      else {
        new_cols <- setdiff(names(expr), colnames(G))
        cols <- c(colnames(G), new_cols)
        tmp <- rep(NA, length(cols))
        names(tmp) <- cols
        intersection <- intersect(cols, names(expr))
        tmp[intersection] <- expr[intersection]
        # Check whether there are genes that need to be added as new columns 
        # of G
        if (length(new_cols) > 0) {
          tmp2 <- matrix(NA, nrow(G), length(new_cols))
          colnames(tmp2) <- new_cols
          rownames(tmp2) <- rownames(G)
          G <- cbind(G, tmp2)
          tmp[new_cols] <- expr[new_cols]
        }
        tmp <- matrix(tmp, nrow=1)
        colnames(tmp) <- cols
        rownames(tmp) <- virus
        # Add the zscores obtained for the current virus as a new row in the 
        # matrix
        G <- rbind(G, tmp)
      }
    }
  }
  return(G)
}

# Getting short names of viruses
viruses <- rownames(B)

# Selecting viruses with gene expression data
sel_viruses <- dir("~/Dropbox/Projetos/Dados/Viruses expression/Differential expression")
sel_viruses <- strsplit(sel_viruses, "_")
sel_viruses <- unlist(sel_viruses)[seq(2, length(sel_viruses)*2, 2)]
sel_viruses <- unlist(strsplit(sel_viruses, "\\."))[seq(1, length(sel_viruses)*2,2)]
sel_viruses <- unique(sel_viruses)
sel_viruses <- intersect(viruses, sel_viruses)

# threshold for selecting significant up or down regulating
thExpr <- 1e-6
G <- geneExpressionMatrix(sel_viruses)
G_orig <- G
# Getting valid gene symbols
symbols <- unique(as.data.frame(org.Hs.egSYMBOL)[,2])

# Selecting valid gene names
sel_genes <- intersect(symbols, colnames(G))

G <- G[,sel_genes]

# Setting non significant differential expression to zero
value <- qnorm(thExpr/2, lower.tail=FALSE)
no_diff_expr <- which(abs(G) <= value)
G[no_diff_expr] <- 0

# Selecting missing values
nas <- which(is.na(G))

# Setting missing values to zero
if (length(nas) > 0)
  G[nas] <- 0

# Removing columns with only zeros
remove_cols <- which(colSums(abs(G), na.rm = TRUE) == 0)
if (length(remove_cols) > 0) {
  G <- G[,-remove_cols]
}

# Making all non-zero values from G equal to 1 (up-regulation) or 
# -1 (down-regulation)
i <- which(G > 0)
G[i] <- 1
i <- which(G < 0)
G[i] <- -1

# Adding viruses with no gene expression data to G
tmp <- matrix(0, length(viruses), ncol(G))
sel_viruses <- viruses
rownames(tmp) <- sel_viruses
colnames(tmp) <- colnames(G)
tmp[intersect(rownames(G), sel_viruses), colnames(G)] <- G[intersect(rownames(G), sel_viruses),]
G <- tmp


