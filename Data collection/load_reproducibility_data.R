# Load 200 embeddings and concatenate the top 100 solutions with lowest cost.
loadConcateMatrices <- function(folder, top=100) {
  files <- dir(folder)
  
  Vr <- NULL
  Dr <- NULL
  Pr <- NULL
  Ur <- NULL
  
  costs <- c()
  for (file in files) {
    load(paste0(folder, file))
    costs[file] <- results$cost[length(results$cost)]
    k <- ncol(results$V)
    Vr <- cbind(Vr, results$V)
    Dr <- cbind(Dr, results$D)
    Pr <- rbind(Pr, results$P)
    Ur <- rbind(Ur, results$U)
  }
  
  if (!is.null(top)) {
    sorted  <<- sort(costs)
    sel <- names(sorted)[1:top]
    is <- which(files %in% sel)
    print(sel)
    sel_pos <- c()
    for (i in is) {
      init <- (i-1)*k + 1
      end <- i*k
      sel_pos <- c(sel_pos, init:end)
    }
    Vr <- Vr[,sel_pos]
    Dr <- Dr[, sel_pos]
    Pr <- Pr[sel_pos, ]
    Ur <- Ur[sel_pos, ]
  }
  
  Vr <<- Vr
  Dr <<- Dr
  Pr <<- Pr
  Ur <<- Ur
  k <<- k
}
