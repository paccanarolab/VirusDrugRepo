# Combining network medicine and collaborative filtering for drug repurposing against viral diseases
This repository contains data and code for running experiments from the manuscript "Combining network medicine and collaborative filtering for drug repurposing against viral diseases", by Suzana de Siqueira Santos, Haixuan Yang, Aldo Galeano, and Alberto Paccanaro

The following files contains the instructions for running the experiments:

- **run_model.R**. R code for running the proposed method, which obtains efficacy predictions for different viruses and drugs.
- **run_evaluation.R**. R code for calculating evaluation metrics for the proposed method and competitors.
- **run_effect_expr_and_adding_viruses.R**. R code for comparing the performance of the method in two scenarios: (1) with gene expression data *versus* without gene expression data, and (2) with all 143 viruses *versus* with only 55 viruses available in the assessment set.
- **run_biological_interpretation.R**. R code for running the analyses of biological interpretation of the learned signatures.
- **run_reproducibility_analysis.R**. R code for running the reproducibility analysis.
- **run_lambda_sensitivity.R**, R code for evaluating the variance of the performance for different values of the hyperparemeters lambda1, lambda2, lambda3, lambda4.

The repository also has two folders, which contain auxiliary code and data for running the experiments:

- **Data collection**. Data and R code for parsing data from the orginal sources.
- **Data analysis**. Processed data and R code with model implementation and auxiliary functions for the remaining experiments (reproducibility analysis, biological interpretation, etc).

Each folder contains a README file describing the files within the folder.

## Dependencies

No additional package is necessary for learning our model. However, if you want to parse the original data and carry out further analyses, you might need to install the following R packages: 

- ROCR
- Polychrome
- ggplot2
- reshape2
- Rtsne
- scatterplot3d

**IMPORTANT**: For dowloading the complete dataset required for biological interpretation and for parsing the orginal data, please install **git lfs**, and initialize it with the command **git lfs install** before cloning the repository.