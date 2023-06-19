
# Auxiliary functions ----------------------------------------------------------

# It computes the cost function.
# nas indicates the missing values
# Here we assume that G has 0 where values are missing.
objectiveFunction <- function(A, B, G, D, V, P, U, lambda1, lambda2, lambda3, 
                              lambda4, alpha=1, nas=NULL) {
    VU <- V%*%U
    if (!is.null(nas))
      VU[nas] <- 0
    if (alpha == 0)
      U <- 0
    return(sum((A-D%*%P)^2) + sum((B-V%*%P)^2) + sum(alpha*(G-VU)^2) + 
             sum(lambda1*(D^2)) + sum(lambda2*(V^2)) + sum(lambda3*(P^2)) + 
             sum(lambda4*(U^2)))
}


# It creates a matrix with n rows and m cols with random values ranging from
# min to max.
# Optionally, it sets the column names or row names.
initializeEmbedding <- function(n, m, min, max, row_names=NULL,
                                col_names=NULL) {
  X <- matrix(NA, n, m)
  X[,] <- runif(length(X), min, max)
  if (!is.null(row_names))
    rownames(X) <- row_names
  if (!is.null(col_names))
    colnames(X) <- col_names
  return(X)
}

# It checks whether the cost function converges to a fix value.
checkStoppingCondition <- function(value, i, eps, max_it) {
  if (!is.null(max_it))
    return(i <= max_it)
  return(value > eps)
}

# Running model ----------------------------------------------------------------

# It learns drug, virus, protein and gene embedding by minimizing the cost
# function with the Adam gradient descent algorithm.
# INPUT:
# - A is a matrix of 0s and 1s indicating which proteins (columns) are perturbed 
# by each drug (rows).
# - B indicates which proteins (columns) are perturbed by each virus (rows).
# - G is a matrix of -1s, 0s and 1s that indicates which genes (columns) are up 
# or down regulated for each virus (row).
# - k is the dimension of the embedding
# - lambda1, lambda2, lamnda3, and lambda4 are parameters for controlling the
# regularization of D, V, P, U, respectively.
# - alpha is a parameter controlling the contribution of the gene expression data
# - eps is used for checking the stopping criterion when max_it is NULL
# - max_it fixes the number of iterations, if it is provided
# - eta, beta1, and beta2 are parameter used by Adam gradient descent algorithm
# - epsilon is small number added for avoiding divisions by 0
# - init_value is a positive number used for defining the range for the initial
# values in the embedding. D, V, and P are initialized with random values from
# a untiform distribution between 0 and init_value. U is initialized with random
# values between -init_value and +init_value.
# OUTPUT
# List containing:
# - scores is a matrix with the final prediction scores
# - D - drug representation in a k-dimension space, with latent features on the
# columns
# - V - virus representation in a k-dimension space, with latent features on the
# columns
# - P - protein representation in a k-dimension space, with latent features on
# the rows
# - U - gene representation in a k-dimension space, with latent features on the
# rows
# - cost is a vector of costs per iteration
gradientDescentAdam <- function(A, B, G, k=15, lambda1=0.1, lambda2=0.1,
                                lambda3=0.1, lambda4=0.1, alpha=1, eps=10e-4, 
                                max_it=NULL, eta=0.001, beta1=0.9, beta2=0.999, 
                                epsilon=1e-8, init_value=0.001) {
  
  # Initializing embeddings with random values
  D <- initializeEmbedding(nrow(A), k, 0, init_value, row_names = rownames(A))
  V <- initializeEmbedding(nrow(B), k, 0, init_value, row_names = rownames(B))
  P <-initializeEmbedding(k, ncol(A), 0, init_value, col_names = colnames(A))
  U <- initializeEmbedding(k, ncol(G), -init_value, init_value,
                           col_names = colnames(G))
  
  # Setting missing values to 0
  nas <- which(is.na(G))
  G[nas] <- 0
  
  
  # Number of iterations
  i <- 1
  
  # Vector of cost per iteration
  cost <- c()
  
  cost[i] <- objectiveFunction(A, B, G, D, V, P, U, lambda1, lambda2, lambda3, 
                               lambda4, alpha=alpha, nas=nas)
  
  # Change of cost function between two iteration
  change <- 1
  
  # Variables used for Adam optimization
  mD <- mV <- mP <- mU <- 0
  vD <- vV <- vP <- vU <- 0

  # Iterations for updating embeddngs until convergence criterion is satisfied  
  while (checkStoppingCondition(change, i,eps,max_it)) {
    print(i)
    # Derivative for updating D
    gradD <- -A%*%t(P) + D%*%P%*%t(P) + lambda1*D
    mD <- beta1*mD + (1-beta1)*gradD
    vD <- beta2*vD + (1-beta2)*(gradD^2)
    mD_ <- mD/(1-beta1^i)
    vD_ <- vD/(1-beta2^i)
    # Updating D
    D <- D - eta*mD_/(sqrt(vD_) + epsilon)
    # Forcing D to be nonnegative
    neg <- which(D < 0)
    D[neg] <- 0
    VU <- V%*%U
    # Adding a mask for ignoring missing values
    if (!is.null(nas))
      VU[nas] <- 0
    # Derivative for updating V
    gradV <- -B%*%t(P) + V%*%P%*%t(P) - alpha*G%*%t(U) + alpha*VU%*%t(U) + lambda2*V
    mV <- beta1*mV + (1-beta1)*gradV
    vV <- beta2*vV + (1-beta2)*(gradV^2)
    mV_ <- mV/(1-beta1^i)
    vV_ <- vV/(1-beta2^i)
    # Updating V
    V <- V - eta*mV_/(sqrt(vV_) + epsilon)
    # Forcing V to be nonnegative
    neg <- which(V < 0)
    V[neg] <- 0
    # Derivative for updating P
    gradP <- -t(D)%*%A + t(D)%*%D%*%P - t(V)%*%B + t(V)%*%V%*%P + lambda3*P
    mP <- beta1*mP + (1-beta1)*gradP
    vP <- beta2*vP + (1-beta2)*(gradP^2)
    mP_ <- mP/(1-beta1^i)
    vP_ <- vP/(1-beta2^i)
    # Updating P
    P <- P - eta*mP_/(sqrt(vP_) + epsilon)
    # Forcing P to be nonnegative
    neg <- which(P < 0)
    P[neg] <- 0
    VU <- V%*%U
    # Addding mask for ignoring missing values
    if (!is.null(nas))
      VU[nas] <- 0
    # # Derivative for updating U
    gradU <- -alpha*t(V)%*%G + alpha*t(V)%*%VU+ lambda4*U
    mU <- beta1*mU + (1-beta1)*gradU
    vU <- beta2*vU + (1-beta2)*(gradU^2)
    mU_ <- mU/(1-beta1^i)
    vU_ <- vU/(1-beta2^i)
    # Updating U
    U <- U - eta*mU_/(sqrt(vU_) + epsilon)
    i <- i + 1
    # Updating cost function
    cost[i] <- objectiveFunction(A, B, G, D, V, P, U, lambda1, lambda2, lambda3, 
                                 lambda4, alpha=alpha, nas=nas)
    # Calculating difference in the cost function
    change <- cost[i-1]-cost[i]
    
  }
  
  scores <- D%*%t(V)
  return(list("scores"=scores, "D"=D, "V"=V, "P"=P, "U"=U, "cost"=cost))
}

