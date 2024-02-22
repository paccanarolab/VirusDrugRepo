# Cosine similarity between rows of matrix A
cosSimilarity <- function(A, norm=TRUE) {
  norms <- sqrt(rowSums(A^2))
  N <- norms%*%t(norms)
  if (!norm)
    N <- 1
  return((A%*%t(A))/N)
}

# Cosine similarity between rows of matrix A and rows of matrix B
cosSimilarityAB <- function(A, B, norm=TRUE) {
  norms1 <- sqrt(rowSums(A^2))
  norms2 <- sqrt(rowSums(B^2))
  
  N <- norms1%*%t(norms2)
  if (!norm)
    N <- 1
  return((A%*%t(B))/N)
}