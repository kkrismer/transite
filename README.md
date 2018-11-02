# Transite

[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)  [![DOI](https://img.shields.io/badge/DOI-10.1101%2F416743-blue.svg)](https://doi.org/10.1101/416743) [![BioC](https://img.shields.io/badge/BioC-1.0.0-brightgreen.svg)](https://doi.org/doi:10.18129/B9.bioc.transite) [![platforms](https://bioconductor.org/shields/availability/3.8/transite.svg)](https://bioconductor.org/packages/release/bioc/html/transite.html#archives)

RNA-binding protein motif analysis

https://transite.mit.edu

Transite is a computational method that allows comprehensive analysis of the regulatory role of RNA-binding proteins in various cellular processes by leveraging preexisting gene expression data and current knowledge of binding preferences of RNA-binding proteins.

Transite and the Transite website were developed in R and C++, using *devtools* to streamline package development, *roxygen2* for documentation, *ggplot2* for visualization, *Rcpp* for C++ integration, *Shiny* as a web framework, *rmarkdown* to generate analysis reports, and *dplyr* for data wrangling.

## About

Transite is a tool to investigate global changes of RNA-binding protein targets.

Transite, a novel computational method that allows cost-effective, time-effective and comprehensive analysis of the regulatory role of RBPs in various cellular processes by leveraging a wealth of preexisting gene expression data and current knowledge of RBP binding preferences. To gain insights into vastly complex processes including the DNA damage response or the immune response, the preliminary step is to calculate the change of mRNA expression levels after stimulus, i.e., the administration of DNA-damaging agents or an immune stimulus, respectively. Based on these results, Transite provides two approaches to investigate inferred mRNA stability changes due to differences in transcript abundance, Transcript Set Motif Analysis and Spectrum Motif Analysis. The former focuses on significantly upregulated and downregulated sets of transcripts and identifies RBPs whose binding sites are overrepresented among those transcripts, whereas the latter approach examines the distribution of RBP binding sites across the entire spectrum of transcripts, sorted according to their fold change, signal-to-noise ratio or any other meaningful measure of differential expression.

## Build status

| Platform | Status |
|------|------|
| Travis CI | [![Travis build status](https://travis-ci.org/kkrismer/transite.svg?branch=master)](https://travis-ci.org/kkrismer/transite) |
| Bioconductor 3.8 (release) | [![BioC release](https://bioconductor.org/shields/build/release/bioc/transite.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/transite/) |
| Bioconductor 3.9 (devel) | [![BioC devel](https://bioconductor.org/shields/build/devel/bioc/transite.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/transite/) |

## Citation

If you use Transite in your research (R package or website), please cite:

**Transite: A computational motif-based analysis platform that identifies RNA-binding proteins modulating changes in gene expression**  
Konstantin Krismer, Shohreh Varmeh, Molly A. Bird, Anna Gattinger, Yi Wen Kong, Thomas Bernwinkler, Daniel A. Anderson, Andreas Heinzel, Brian A. Joughin, Ian G. Cannell, and Michael B. Yaffe  
bioRxiv 416743; doi: https://doi.org/10.1101/416743

## Funding

The development of this package was supported by scholarships of the Marshall Plan Foundation and the Austrian Federal Ministry for Education, National Institutes of Health (NIH) grants R01-ES015339, R35-ES028374, U54-CA112967, the Charles and Marjorie Holloway Foundation, and a Starr Cancer Consortium Award I9-A9-077.
