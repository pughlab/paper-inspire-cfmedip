# INSPIRE cfMeDIP-seq Paper Reproducible Report

In the interest of transparency, openness, and reproducibility, this report contains the analysis code which generates all of the figures and numerical results for the manuscript entitled *Early changes in tumor-agnostic cell-free methylomes and fragmentomes predict outcomes in patients with solid tumors treated with pembrolizumab*. In it, we report on the results of cell-free methylated DNA immunoprecipitation and sequencing (cfMeDIP-seq) from the investigator-initiated phase 2 study of pembrolizumab immunological response evaluation (INSPIRE).

Using Rstudio or the knitr package in R, all results reported in the manuscript can be regenerated in full using this RMarkdown file. All of the processed, deidentified data are stored in `figure_data.Rds` as dataframes and are attached to the environment in the setup block. Our nomenclature is variable names starting with `df` to denote source data tables and `plotdata` for processed dataframes used as intermediates to generate figures. Variables that start with `val` denote numeric values, many of which are reported in the final results.

