load("Data analysis/Input/input_matrices.RData")

source("Data analysis/model.R")

# Runs the model. The output of the function is a list containing the predicition scores,
# and the signatures of drugs (D), viruses (V), proteins (P), and genes (U).
results <- gradientDescentAdam(A, B, G, max_it=30)

# Save prediction scores
write.csv(results$scores, "predictions.csv", quote=FALSE, row.names=TRUE)