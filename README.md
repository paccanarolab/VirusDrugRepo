# Combining network medicine and collaborative filtering for drug repurposing against viral diseases
This repository contains data and code for running experiments from the manuscript "Combining network medicine and collaborative filtering for drug repurposing against viral diseases", by Suzana de Siqueira Santos, Haixuan Yang, Aldo Galeano, and Alberto Paccanaro

It contains 6 files:

- **run_model.R**. R code for running the proposed method, which obtains efficacy predictions for different viruses and drugs.
- **run_evaluation.R**. R code for calculating evaluation metrics for the proposed method and competitors.
- **run_effect_expr_and_adding_viruses.R**. R code for comparing the performance of the method in two scenarios: (1) with gene expression data *versus* without gene expression data, and (2) with all 143 viruses *versus* with only 55 viruses available in the assessment set.
- **run_biological_interpretation.R**. R code for running the analyses of biological interpretation of the learned signatures.
- **run_reproducibility_analysis.R**. R code for running the reproducibility analysis.
- **run_lambda_sensitivity.R**, R code for evaluating the variance of the performance for different values of the hyperparemeters lambda1, lambda2, lambda3, lambda4.

And two folders: 

- **Data collection**. Use the R code provided in this folder if you want to parse data from the orginal sources.
- **Data analysis**. The R code provided in this folder constains the model implementation and further experiments

Each folder contains a README file describing the files within the folder.

## Dependencies

No additional package is necessary for running the proposed method. However, if you want to parse the original data and perform further analyses, you might need to install the following R packages: 

- ROCR
- Polychrome
- ggplot2
- reshape2
- Rtsne
- scatterplot3d
