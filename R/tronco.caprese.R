#### tronco.caprese.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export tronco.caprese
#' @title runs CAPRESE algorithm
#'
#' @description
#' \code{tronco.caprese} run CAPRESE algorithm on the specified dataset. 
#'
#' @details
#' \code{tronco.caprese} infer the topology and compute the confidence measures defined in \code{confidence}.
#' 
#' @param dataset The input dataset. Type: dataframe. The dataset given as input is loaded by the \code{data} function.
#' @seealso \code{\link{data}}
#' @param lambda the real positive value of the shrinkage coefficient, required to range in [0, 1]. Its default value is 0.5.
#' @param do.estimation run CAPRESE algorithm and estimates the error rates. Type: boolean, dafault: TRUE.
#' @param verbose run CAPRESE algorithm with verbose output. Type: boolean, dafault: FALSE.
#' @return an object containing the reconstructed topology and confidence values.
#' 
tronco.caprese <- function(dataset = NA, lambda = 0.5, do.estimation = TRUE, verbose = FALSE) {
	if(is.na(dataset)) {
		dataset = settings$data.values$dataset;
	}
	if(lambda < 0 || lambda > 1) {
		stop("The value of the lambda shrinkage parameter has to be in [0:1]!",call.=FALSE);
	}
	#reconstruct the topology with CAPRESE
	topology = caprese.fit(dataset,lambda,do.estimation,verbose);
	#set the labels to the reconstructions
	#labels for the dataset
    colnames(topology$dataset) <- settings$visualization$labels;
    #labels for the reconstruction
    rownames(topology$adj.matrix) <- settings$visualization$labels;
    colnames(topology$adj.matrix) <- settings$visualization$labels;
    rownames(topology$parents.pos) <- settings$visualization$labels;
    rownames(topology$confidence) <- settings$visualization$labels;
    colnames(topology$confidence) <- settings$visualization$labels;
    #labels for the probabilities
    rownames(topology$probabilities$marginal.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$joint.probs) <- settings$visualization$labels;
    colnames(topology$probabilities$joint.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$conditional.probs) <- settings$visualization$labels;
    if(do.estimation==TRUE) {
		rownames(topology$probabilities$estimated.marginal.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$estimated.joint.probs) <- settings$visualization$labels;
		colnames(topology$probabilities$estimated.joint.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$estimated.conditional.probs) <- settings$visualization$labels;
	}
	else {
		topology$probabilities$estimated.marginal.probs = matrix();
		topology$probabilities$estimated.joint.probs = matrix();
		topology$probabilities$estimated.conditional.probs = matrix();
		topology$error.rates = list(error.fn=numeric(),error.fp=numeric());
	}
    #create the new object topology
    topology.obj <- new("topology",
		#the input valid dataset
		dataset = topology$dataset,
		#the observed and estimated probabilities
		marginal.probs = topology$probabilities$marginal.probs,
		joint.probs = topology$probabilities$joint.probs,
		conditional.probs = topology$probabilities$conditional.probs,
		estimated.marginal.probs = topology$probabilities$estimated.marginal.probs,
		estimated.joint.probs = topology$probabilities$estimated.joint.probs,
		estimated.conditional.probs = topology$probabilities$estimated.conditional.probs,
		#the observed and estimated probabilities for the prima facie topology
		pf.marginal.probs = matrix(),
		pf.joint.probs = matrix(),
		pf.conditional.probs = matrix(),
		pf.estimated.marginal.probs = matrix(),
		pf.estimated.joint.probs = matrix(),
		pf.estimated.conditional.probs = matrix(),
		#the reconstructed topology
		adj.matrix = topology$adj.matrix,
		parents.pos = topology$parents.pos,
		#the reconstructed prima facie topology
		pf.adj.matrix = matrix(),
		pf.parents.pos = matrix(),
		#the estimated error rates for the reconstructed and prima facie topologies
		error.fp = topology$error.rates$error.fp,
		error.fn = topology$error.rates$error.fn,
		pf.error.fp = numeric(),
		pf.error.fn = numeric(),
		#confidence in terms of prima facie scores
		confidence.scores = list(pr.scores=topology$confidence),
		#confidence in terms of bootstrap
		confidence.np = list(),
		pf.confidence.np = list(),
		edge.confidence.np = matrix(),
		pf.edge.confidence.np = matrix(),
		bootstrap.settings.np = list(),
		bootstrap.np = FALSE,
		confidence.p = list(),
		pf.confidence.p = list(),
		edge.confidence.p = matrix(),
		pf.edge.confidence.p = matrix(),
		bootstrap.settings.p = list(),
		bootstrap.p = FALSE,
		#parameters
		parameters = topology$parameters,
		algorithm = "CAPRESE")
	#the reconstruction has been completed
	cat(paste("Executed CAPRESE algorithm with shrinkage coefficient:",lambda,"\n"));
	if(do.estimation==TRUE) {
		cat(" Estimated false positives error rate: ",topology.obj@error.fp);
		cat("\n Estimated false negative error rate: ",topology.obj@error.fn);
	}
    #save the reconstruction to the global workspace and return it
    assign("topology",topology.obj,envir=.GlobalEnv);
    return(topology.obj);
}

#### end of file -- tronco.caprese.R
