# Transite

[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)  [![DOI](https://img.shields.io/badge/DOI-10.1101%2F416743-blue.svg)](https://doi.org/10.1101/416743) [![BioC](https://img.shields.io/badge/BioC-1.6.2-brightgreen.svg)](https://doi.org/doi:10.18129/B9.bioc.transite) [![platforms](https://bioconductor.org/shields/availability/3.11/transite.svg)](https://bioconductor.org/packages/release/bioc/html/transite.html#archives) [![Coverage Status](https://coveralls.io/repos/github/kkrismer/transite/badge.svg?branch=master)](https://coveralls.io/github/kkrismer/transite?branch=master)

RNA-binding protein motif analysis

https://transite.mit.edu

Transite is a computational method that allows comprehensive analysis of the regulatory role of RNA-binding proteins in various cellular processes by leveraging preexisting gene expression data and current knowledge of binding preferences of RNA-binding proteins.

Transite and the Transite website were developed in R and C++, using *devtools* to streamline package development, *Travis CI* as continuous integration service, *Git* for version control, *roxygen2* for documentation, *ggplot2* for visualization, *Rcpp* for C++ integration, *Shiny* as web framework, *rmarkdown* and *knitr* to generate analysis reports, and *dplyr* for data wrangling.

## About

Transite is a tool to investigate global changes of RNA-binding protein targets.

Transite, a novel computational method that allows cost-effective, time-effective and comprehensive analysis of the regulatory role of RBPs in various cellular processes by leveraging a wealth of preexisting gene expression data and current knowledge of RBP binding preferences. To gain insights into vastly complex processes including the DNA damage response or the immune response, the preliminary step is to calculate the change of mRNA expression levels after stimulus, i.e., the administration of DNA-damaging agents or an immune stimulus, respectively. Based on these results, Transite provides two approaches to investigate inferred mRNA stability changes due to differences in transcript abundance, Transcript Set Motif Analysis and Spectrum Motif Analysis. The former focuses on significantly upregulated and downregulated sets of transcripts and identifies RBPs whose binding sites are overrepresented among those transcripts, whereas the latter approach examines the distribution of RBP binding sites across the entire spectrum of transcripts, sorted according to their fold change, signal-to-noise ratio or any other meaningful measure of differential expression.

## Installation

The *transite* package is part of Bioconductor since release 3.8. To install it on your system, enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("transite")
```

The above method is recommended for most users, as it is the most recent stable version.

The development version can be installed directly from this repository (most recent, least stable):

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("kkrismer/transite")
```

Alternatively, the development version can be installed using the development branch of Bioconductor:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(version = "devel"")
BiocManager::install("transite")
```
Note: For most use cases it is not necessary to install the `transite` package locally, as a substantial part of its functionality is offered as an online service at https://transite.mit.edu.

## Build status

| Platform | Status |
|------|------|
| Travis CI | [![Travis build status](https://travis-ci.org/kkrismer/transite.svg?branch=master)](https://travis-ci.org/kkrismer/transite) |
| Bioconductor 3.11 (release) | [![BioC release](https://bioconductor.org/shields/build/release/bioc/transite.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/transite/) |
| Bioconductor 3.12 (devel) | [![BioC devel](https://bioconductor.org/shields/build/devel/bioc/transite.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/transite/) |

## Citation

If you use Transite in your research (R package or website), please cite:

**Transite: A computational motif-based analysis platform that identifies RNA-binding proteins modulating changes in gene expression**  
Konstantin Krismer, Molly A. Bird, Shohreh Varmeh, Erika D. Handly, Anna Gattinger, Thomas Bernwinkler, Daniel A. Anderson, Andreas Heinzel, Brian A. Joughin, Yi Wen Kong, Ian G. Cannell, and Michael B. Yaffe  
Cell Reports, Volume 32, Issue 8, 108064; DOI: https://doi.org/10.1016/j.celrep.2020.108064

## Funding

The development of this package was supported by scholarships of the Marshall Plan Foundation and the Austrian Federal Ministry for Education, National Institutes of Health (NIH) grants R01-ES015339, R35-ES028374, U54-CA112967, the Charles and Marjorie Holloway Foundation, the MIT Center for Precision Cancer Medicine, and a Starr Cancer Consortium Award I9-A9-077.
