TRONCO (TRanslational ONCOlogy)
===============================

| Branch              | Stato CI      |
|---------------------|---------------|
| master | [![Build Status](https://travis-ci.org/BIMIB-DISCo/TRONCO.svg?branch=master)](https://travis-ci.org/BIMIB-DISCo/TRONCO) |
| development | [![Build Status](https://travis-ci.org/BIMIB-DISCo/TRONCO.svg?branch=development)](https://travis-ci.org/BIMIB-DISCo/TRONCO) |


The **TRONCO** (*TR*anslational *ONCO*logy) **R** package collects algorithms to infer *progression models* via the approach of Suppes-Bayes Causal Network, both from an ensemble of tumors (cross-sectional samples) and within an individual patient (multi-region or single-cell samples). 

The package provides parallel implementation of algorithms that process binary matrices where each row represents a tumor sample and each column a single-nucleotide or a structural variant driving the progression; a 0/1 value models the absence/presence of that alteration in the sample. 

The tool can import data from plain, *MAF* or *GISTIC* format files, and can fetch it from the cBioPortal for cancer genomics. Functions for  data manipulation and visualization are provided, as well as functions to import/export such data to other bioinformatics  tools for, e.g, clustering or detection of mutually exclusive alterations. Inferred models can be visualized and tested for their confidence via bootstrap and cross-validation. 

TRONCO is used for the implementation of the Pipeline for Cancer Inference (PICNIC). 
