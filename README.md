# INSPIRE cfMeDIP-seq Paper Reproducible Report

In the interest of transparency, openness, and reproducibility, this report contains the analysis code which generates all of the figures and numerical results for the manuscript entitled *Early changes in tumor-agnostic cell-free methylomes and fragmentomes predict outcomes in patients with solid tumors treated with pembrolizumab*. In it, we report on the results of cell-free methylated DNA immunoprecipitation and sequencing (cfMeDIP-seq) from the investigator-initiated phase 2 study of pembrolizumab immunological response evaluation (INSPIRE).

## Instructions

Using Rstudio or the knitr package in R, all results reported in the manuscript can be regenerated in full using the RMarkdown file, `INSPIRE-cfMeDIP-Paper-Code.Rmd`.

All of the processed, deidentified data are stored in `figure_data.Rds` as dataframes and are attached to the environment in the setup block. Our nomenclature is variable names starting with `df` to denote source data tables and `plotdata` for processed dataframes used as intermediates to generate figures. Variables that start with `val` denote numeric values, many of which are reported in the final results.

To view the results of the analysis, download `INSPIRE-cfMeDIP-Paper-Code.html`.

To run the analysis yourself, follow these steps:

1. Install R - https://www.r-project.org/
2. Install Rstudio - https://posit.co/products/open-source/rstudio/
3. Clone this repository in the terminal with `git clone https://github.com/pughlab/paper-inspire-cfmedip.git` **OR** download this repository as a ZIP file.
4. Open `INSPIRE-cfMeDIP-Paper-Code.Rmd` in Rstudio.
5. Install `knitr`, `rmdformats` and all of the dependency packages which are called by `library()` function in the `setup` block.
6. Run all blocks or knit the R Markdown document.
