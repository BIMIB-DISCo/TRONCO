[![Github Issues](https://img.shields.io/github/issues/badges/shields.svg?style=flat)](https://github.com/BIMIB-DISCo/TRONCO/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/zachmayer/caretEnsemble/pulls.svg?style=flat)](https://github.com/BIMIB-DISCo/TRONCO/issues/pulls)

TRONCO (TRanslational ONCOlogy)
===============================

**TRONCO** is a **R** package which collects algorithms to infer *progression models* from Bernoulli 0/1 profiles of genomic alterations across a tumor sample. 

Such profiles are usually visualised as a binary input matrix where each row represents a patient’s sample (e.g., the result of a sequenced tumor biopsy), and each column an event relevant to the progression (a certain type of somatic mutation, a focal or higher-level chromosomal copy number alteration etc.); a 0/1 value models the absence/presence of that alteration in the sample. 

In this version of **TRONCO** such profiles can be readily imported by boolean matrices and *MAF* or *GISTIC* files. The package provides various functions to editing, visualise and subset such data, as well as functions to query the Cbio portal for cancer genomics. 

This version of TRONCO comes with the parallel implementations the **CAPRESE**  [*PLoS ONE 9(12): e115570*] and **CAPRI** [*Bioinformatics, doi:10.1093/bioinformatics/btv296*] algorithms to infer possible progression models arranged as trees, or general direct acyclic graphs. Bootstrap functions to assess the parametric, non-prametric and statistical confidence of every inferred model are also provided. The package comes with some data available as well, which include the dataset of *Atypical Chronic Myeloid Leukemia samples* provided by Piazza et al., Nat. Genet., 45 (2013), and examples.
