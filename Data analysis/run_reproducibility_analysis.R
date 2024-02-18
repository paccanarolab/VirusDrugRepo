source("Data analysis/reproducibility.R")
source("Data collection/load_reproducibility_data.R")

loadConcateMatrices("Data analysis/Embedding/Different initializations - alpha = 1/")

repV <- repD <- repP <- repU <- list()

set.seed(0)
n_kmeans <- 100 # number of times kmeans is executed 
seeds <- runif(n_kmeans, 0, n_kmeans)

# Run clustering algorithm 100 times
for (i in 1:n_kmeans) {
  repV[[i]] <- repAnalysis(Vr, k, seed=seeds[i])
  repD[[i]] <- repAnalysis(Dr, k, seed=seeds[i])
  repP[[i]] <- repAnalysis(t(Pr), k, seed=seeds[i])
  repU[[i]] <- repAnalysis(t(Ur), k, seed=seeds[i])
  
}

# Compute the mean silhouette for each run of the clustering algorithm
mean_silhoutteD <- meanSilhouette(repD)
mean_silhoutteV <- meanSilhouette(repV)
mean_silhoutteP <- meanSilhouette(repP)
mean_silhoutteU <- meanSilhouette(repU)

# Select the best run of the clustering algorithm
bestV <- repV[[which.max(mean_silhoutteV)]]
bestD <- repD[[which.max(mean_silhoutteD)]]
bestP <- repP[[which.max(mean_silhoutteP)]]
bestU <- repU[[which.max(mean_silhoutteU)]]

# Plot the silhouette obtained for the virus embedding
pdf("Data analysis/Plots/silhouttePlotV.pdf")
plotSilhouette(bestV, computeSil=FALSE)
dev.off()

# Plot the silhouette obtained for the drug embedding
pdf("Data analysis/Plots/silhouttePlotD.pdf")
plotSilhouette(bestD, computeSil=FALSE)
dev.off()

# Plot the silhouette obtained for the virus embedding
pdf("Data analysis/Plots/silhouttePlotP.pdf")
plotSilhouette(bestP, computeSil=FALSE)
dev.off()

# Plot the silhoette obtained for the gene embedding
pdf("Data analysis/Plots/silhouttePlotU.pdf")
plotSilhouette(bestU, computeSil=FALSE)
dev.off()
