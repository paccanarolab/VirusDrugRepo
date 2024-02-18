require(ROCR)

source("Data collection/load_evaluation_data.R")

computeAucAprLabels <- function (scores, labels) {
  
  scores <- matrix(scores, ncol=1)
  auc = NA
  apr = NA
  try({x   <- prediction(scores, labels)})
  try({ROC <- performance(x, "tpr", "fpr")})
  try({PRC <- performance(x, "prec", "rec")})
  
  # for (m in 1:ncol(drug_scores)) {
  try({f <- approxfun(data.frame(ROC@x.values[[1]], ROC@y.values[[1]]))}) 
  x <- seq(0, 1, 1/10000)
  try({y <- unlist(lapply(x, f))})
  try({auc <- trapezoidSum(x,y)})
  #aucs[m] <- integrate(f, 0, 1)$value
  #}
  try({f <- approxfun(data.frame(PRC@x.values[[1]], PRC@y.values[[1]]))}) 
  x <- seq(0, 1, 1/10000)
  try({y <- unlist(lapply(x, f))})
  try({apr <- trapezoidSum(x,y)})
  
  return(list("auc"=auc, "apr"=apr))
}


evaluate <- function(scores, viruses, sel_drugs=NULL) {
  drugs <- rownames(scores)
  if (!is.null(sel_drugs))
    drugs <- sel_drugs
  loadDrugVirusAssociations()
  aucs <- c()
  aprs <- c()
  topns <- c(1, 3, 5, 10, 20, 30, 40, 50, 100, 150, 200, 500, 1000, 1500, 2000)
  recall_per_virus <- matrix(0, length(viruses), length(topns))
  rownames(recall_per_virus) <- viruses
  colnames(recall_per_virus) <- as.character(topns)
  
  recall <- rep(0, length(topns))
  names(recall) <- topns
  for (virus in viruses) {
    labels <- drug_virus[,virus]
    sc <- rep(0, length(labels))
    names(sc) <- rownames(drug_virus)
    sc[rownames(scores)] <- scores[,virus]
    if (sum(is.na(sc)) > 0) {
      aucs[virus] <- NA
      aprs[virus] <- NA
    }
    else {
      auc_apr <- computeAucAprLabels(sc, matrix(labels, ncol=1))
      aucs[virus] <- auc_apr$auc
      aprs[virus] <- auc_apr$apr
      sorted_scores <- sort(scores[, virus], decreasing=TRUE)
      for (top in topns) {
        recall_per_virus[virus, as.character(top)] <- length(which(names(sorted_scores[1:top]) %in% sel_drugs))
        recall[as.character(top)] <- recall[as.character(top)] + recall_per_virus[virus, as.character(top)]
        recall_per_virus[virus, as.character(top)] <- recall_per_virus[virus, as.character(top)]/length(sel_drugs)
      }
    }
  }
  total_associations <- sum(drug_virus)
  recall <- recall/sum(drug_virus)
  output <- list("auc" = aucs, "apr"=aprs, "recall"=recall, "recall_per_virus"=recall_per_virus)
  return(output)
}

evaluateMethods <- function(scores_list, viruses, sel_drugs) {
  methods <- names(scores_list)
  aucs <- matrix(NA, length(viruses), length(methods))
  rownames(aucs) <- viruses
  colnames(aucs) <- methods
  apr <- aucs 
  recall_per_virus_150 <- aucs
  recall_per_virus_20 <- aucs
  recall_per_virus_50 <- aucs
  recall_per_virus_100 <- aucs
  recall_per_virus_200 <- aucs
  for (method in methods) {
    sc <- scores_list[[method]][sel_drugs, viruses]
    if (length(viruses) == 1) {
      sc <- matrix(sc, ncol=1)
      rownames(sc) <- sel_drugs
      colnames(sc) <- viruses
    }
    eval <- evaluate(sc, viruses, sel_drugs = sel_drugs)
    aucs[viruses, method] <- eval$auc
    apr[viruses, method] <- eval$apr
    recall_per_virus_150[viruses, method] <- eval$recall_per_virus[, "150"]
    recall_per_virus_20[viruses, method] <- eval$recall_per_virus[, "20"]
    recall_per_virus_50[viruses, method] <- eval$recall_per_virus[, "50"]
    recall_per_virus_100[viruses, method] <- eval$recall_per_virus[, "100"]
    recall_per_virus_200[viruses, method] <- eval$recall_per_virus[, "200"]
    
  }
  return(list("auc"=aucs, "apr"=apr, "recall"=recall_per_virus_150, 
              "recall_20"=recall_per_virus_20, "recall_50"=recall_per_virus_50,
              "recall_100"=recall_per_virus_100, "recall_200"=recall_per_virus_200))
  
}

evaluateFiles <- function(dir, files) {
  loadDrugVirusAssociations()
  sel_viruses <- which(colSums(drug_virus) == 0)
  sel_viruses <- colnames(drug_virus)[sel_viruses]
  aucs <- matrix(NA, length(sel_viruses), length(files))
  rownames(aucs) <- sel_viruses
  colnames(aucs) <- files
  
  recall <- aucs
  
  for (file in files) {
    scores <- getPredictions(paste0(dir, file))
    eval <- evaluate(scores, sel_viruses)
    aucs[, file] <- eval$auc[viruses_hvidb_eval]
    recall[, file] <- eval$recall_per_virus[viruses_hvidb_eval, "150"]
  }
  
  return(list("auc"=aucs, "recall"=recall))
}

