This folder contains the files:

    DrugBank_ATC.tsv - text file with drug ATC categories
    drug_chem_sim.RData - Rdata file used for loading variable chem_sim, with the chemical similarity between drugs
    drug_kegg_pathways.RData -  R data file with kegg pathways altered in cell lines treated by drugs
    drug_pathways.RData - Rdata file with drugbank pathways
    kegg_pathways.RData - Rdata with sets of proteins (entrez ids and gene symbols) for each KEGG pathway
    kernel.RData - Rdata with square matrix corresponding to the p-step random walk on the PPI
    pvalues_drug_virus_kegg.RData - pvalue of the Wilcoxon test comparing the cosine similarity between viruses and drugs 
    that share the same pathways and that that do not, for each KEGG pathway altered in infected cell lines or treated cell
    lines.
    drug_targets_Gysi.RData - RData file used for loading variable drug_targets, with the known drug-target associations
    drug_virus.RData - RData file with evaluation data (known drug-virus associations)
    sim_bp_entrez.RData - RData file with (BP - Biological Processes) semantic similarity between drug targets 
    sim_cc_entrez.RData - RData file with (CC - Cellular Component) semantic similarity between drug targets 
    sim_mf_entrez.RData - RData file with (MF - Molecular function) semantic similarity between drug targets 
    virus_diff_expressed_pathways.RData - Rdata file with kegg pathways altered in cell lines infected by viruses
    virus_families.RData - RData with virus familiies
    virus_host_Gysi.RData - RData with virus-host associations from HVIDB.

And the subfolder:
    'Gene expression' - folder with gene expression data for 34 viruses. For each virus it has a file with a
    z-score measuring the differential expression.
    