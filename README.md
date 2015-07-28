TRONCO (TRanslational ONCOlogy)
===============================

**TRONCO** is a **R** package which collects algorithms to infer *progression models* from Bernoulli 0/1 profiles of genomic alterations across a tumor sample. 

Such profiles are usually visualized as a binary input matrix where each row represents a patient’s sample (e.g., the result of a sequenced tumor biopsy), and each column an event relevant to the progression (a certain type of somatic mutation, a focal or higher-level chromosomal copy number alteration, etc.); a 0/1 value models the absence/presence of that alteration in the sample. 

In this version of **TRONCO** such profiles can be readily imported by boolean matrices and *MAF* or *GISTIC* files. The package provides various functions to editing, visualize and subset such data, as well as functions to query the cBioPortal for cancer genomics. 

In the current version, **TRONCO** provides parallel implementations of **CAPRESE**  [*PLoS ONE 9(12): e115570*] and **CAPRI** [*Bioinformatics, doi:10.1093/bioinformatics/btv296*] algorithms to infer progression models arranged as trees or general direct acyclic graphs. Bootstrap procedures to assess the non-prametric and statistical confidence of the inferred models are also provided. The package comes with example data available, which include the dataset of *Atypical Chronic Myeloid Leukemia samples* by Piazza et al., Nat. Genet., 45 (2013). 
