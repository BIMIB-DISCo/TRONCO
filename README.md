TRONCO (TRanslational ONCOlogy)
===============================

[![Actions Status](https://github.com/BIMIB-DISCo/TRONCO/workflows/check-master/badge.svg)](https://github.com/BIMIB-DISCo/TRONCO/actions?query=workflow%3Acheck-master)
[![Actions Status](https://github.com/BIMIB-DISCo/TRONCO/workflows/check-development/badge.svg)](https://github.com/BIMIB-DISCo/TRONCO/actions?query=workflow%3Acheck-development)

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/BIMIB-DISCo/TRONCO/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/BIMIB-DISCo/TRONCO/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/BIMIB-DISCo/TRONCO/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/BIMIB-DISCo/TRONCO/actions/workflows/check-bioc.yml) |


The **TRONCO** (*TR*anslational *ONCO*logy) **R** package collects algorithms to infer *progression models* via the approach of Suppes-Bayes Causal Network, both from an ensemble of tumors (cross-sectional samples) and within an individual patient (multi-region or single-cell samples). 

The package provides parallel implementation of algorithms that process binary matrices where each row represents a tumor sample and each column a single-nucleotide or a structural variant driving the progression; a 0/1 value models the absence/presence of that alteration in the sample. 

The tool can import data from plain, *MAF* or *GISTIC* format files, and can fetch it from the cBioPortal for cancer genomics. Functions for  data manipulation and visualization are provided, as well as functions to import/export such data to other bioinformatics  tools for, e.g, clustering or detection of mutually exclusive alterations. Inferred models can be visualized and tested for their confidence via bootstrap and cross-validation. 

TRONCO is used for the implementation of the Pipeline for Cancer Inference (PICNIC). 
