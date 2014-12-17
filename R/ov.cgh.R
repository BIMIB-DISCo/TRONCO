#### ov.cgh.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @name ov.cgh
#' @title Ovarian cancer CGH data
#' @description 
#' This dataset is obtained using the comparative genomic hybridization technique (CGH) 
#' on samples, i.e., patiens, affected by papillary serous cystadenocarcinoma of the 
#' ovary. Only the seven more frequent events are given.
#' @docType data
#' @usage data.load("CGH")
#' @format
#' A data frame of 87 observations (patiens) on 7 variables (mutations).
#' @source \url{http://www.ncbi.nlm.nih.gov/sky/}
#' @details
#' The CGH technique uses fluorescent staining to detect abnormal (increased or 
#' decreased) number of DNA copies. Often the results are reported as a gain or loss on 
#' a certain arm, without further distinction for specific regions. It is common to 
#' denote a change in DNA copy number on a specific chromosome arm by prefixing a "-" 
#' sign for loss and a "+" for gain. Thus, say, -3q denotes abnormally low DNA copy 
#' number on the q arm of the 3rd chromosome.
NULL

#### end of file -- ov.cgh.R
