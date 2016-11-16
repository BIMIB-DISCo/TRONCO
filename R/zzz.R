#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

#' @useDynLib TRONCO
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
	hypotheses.env <- new.env(parent = emptyenv())
	assign('hypotheses.env', hypotheses.env, asNamespace(pkgname))
}

.onUnload <- function (libpath) {
  library.dynam.unload("TRONCO", libpath)
}

#### end of file -- as.functions.R
