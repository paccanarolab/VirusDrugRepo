This folder contains the file 
    selectedEmbedding.RData - RData file with drugs, virus, protein and gene signatures learned by our method. We selected the signatures with the lowest cost function at convergence
    
and the folders 
'Different initializations - alpha = 0' - embbeddings obtained by our method with 200 different initializations without using gene expression data
'Different initializations - alpha = 1' -  embbeddings obtained by our method with 200 different initializations 
'Different initializations - Assessment set' - embbeddings obtained by our method with 200 different initializations using a smaller dataset (Assessement set)
'Different lambdas' - embeddings obtained by our method using different values of the hyperparameters lambda1, lambda2, lambda3, and lambda4