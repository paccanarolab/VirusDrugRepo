
loadPredictionsCompetitors <- function() {
  load("Data analysis/Competitors/predictions_competitors.RData")
  santos <<- santos
  li_et_al <<- li_et_al
  save_runner <<- save_runner_adj
}

loadDrugVirusAssociations <- function() {
  load("Data collection/Data/drug_virus.RData")
  drug_virus <<- drug_virus
}

loadPredctionsLambdas <- function() {
  folder <- "Data analysis/Embedding/Different lambdas/"
  files <- dir(folder)
  names <- unlist(strsplit(files, "="))
  names <- names[seq(2, length(names), 2)]
  names <- unlist(strsplit(names,".RData"))
  names <- gsub("_", "/", names)
  names(files) <- names
  scores_lambdas_diff <- list()
  i <- grep("/", names)
  for (name in names[i])
    scores_lambdas_diff[[name]] <- getPredictions(paste0(folder, files[name]))
  scores_lambdas <- list()
  for (name in names[-i])
    scores_lambdas[[name]] <- getPredictions(paste0(folder, files[name]))
  i <- order(as.numeric(names(scores_lambdas)))
  scores_lambdas <- scores_lambdas[i]
  scores_lambdas <<- scores_lambdas
  scores_lambdas_diff <<- scores_lambdas_diff
}

getPredictions <- function(file) {
  load(file)
  scores <- results$D%*%t(results$V)
  return(scores)
}