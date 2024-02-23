library("GMKMcharlie")

source("Data analysis/cosine_similarity.R")


silhouettePerCluster <- function(clusters, silhouttes) {
  res <- list()
  for (i in 1:length(clusters)) {
    res[[i]] <- silhouttes[clusters[[i]]$clusterMember]
  }
  return(res)
}

mapClusters <- function(i, clusters) {
  for (j in 1:length(clusters)) {
    if (i %in% clusters[[j]]$clusterMember)
      return(j)
  }
}

mapAllItemsCluster <- function(clusters, nitems) {
  map <- c()
  for (i in 1:nitems) {
    map[i] <- mapClusters(i, clusters)
  }
  return(map)
}


listOfClusters <- function(clusters) {
  elements <- c()
  for (i in 1:length(clusters)) {
    elements[clusters[[i]]$clusterMember] <- i
  }
  return(elements)
}

# Cluster the columns of the concatenated matrix into k clusters
clusterComponents <- function(W, k, init="++", seed=0) {
  # Components are on the columns
  set.seed(seed)
  initial_node <- sample(1:ncol(W), 1)
  centers <- GMKMcharlie::KMppIni(X = W, K = k, seed=seed, stochastic=TRUE, 
                                  firstSelection = initial_node)
  values <- W[, centers]
  clusters <- GMKMcharlie::KM(X = W, centroid = values, minkP = "cosine")
  return(clusters)
}

# Compute the reproducibility score for a given concatenated matrix
repAnalysis <- function(W, k, init="++", seed=0) { 
  cluster <- clusterComponents(W, k, init=init, seed = seed)
  dist <- 1 - cosSimilarity(t(W))
  i <- which(dist < 0)
  if (length(i) > 0)
    dist[i] <- 0
  sil <- cluster::silhouette(listOfClusters(cluster), dmatrix=dist) 
  s <- sil[, 'sil_width']
  f <- function(x){return(length(x$clusterMember))}
  ncluster <- unlist(lapply(cluster, f))
  scluster <- silhouettePerCluster(cluster, s)
  return(list("cluster"=cluster, "dist"=dist, "silhouette"=s, "sil"=sil, 
              "ncluster"=ncluster, "scluster"=scluster))
}

# For each experiment compute the mean silhouette across the clustes 
meanSilhouette <- function(rep) {
  tmp <- c()
  for (i in 1:length(rep))
    tmp[i] <- mean(rep[[i]]$silhouette)
  return(tmp)
}

# Plot the silhouettes
plotSilhouette <- function(rep) {
  sizes <- rep$ncluster
  cols <- c()
  colors <- rep(2:8, rep(length(sizes)))
  for (i in 1:length(sizes)) cols <- c(cols, rep(colors[i], sizes[i]))
  sil <- rep$sil
  plot(sil, col=cols, main="Silhouette plot")
}
