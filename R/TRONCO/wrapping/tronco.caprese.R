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
#' \code{tronco.caprese} run CAPRESE algorithm on the specified inputs. 
#'
#' @details
#' \code{tronco.caprese} infer the best tree topology and compute the confidence measures defined in \code{confidence}.
#' 
#' @param data The input data. Type: dataframe. The dataset provided inside data is assumed to be valid.
#' @seealso \code{\link{data}}
#' @param lambda a real positive value defining the shrinkage coefficient, required to range in [0, 1]. Its default value is 0.5.
#' @param do.estimation run CAPRESE algorithm and estimates the error rates. Type: boolean, dafault: FALSE.
#' @return a list saving the reconstructed tree topology and the confidence values.
#' 
tronco.caprese <- function(data, lambda = 0.5, do.estimation = FALSE) {
	#check for the inputs to be correct
	if(is.null(data) || is.null(data$genotypes)) {
		stop("The dataset given as input is not valid.");
	}
	if(lambda < 0 || lambda > 1) {
		stop("The value of the shrinkage parameter lambda has to be in [0:1]!",call.=FALSE);
	}
	#reconstruct the topology with CAPRESE
	cat(paste("Running CAPRESE algorithm with shrinkage coefficient:",lambda,"\n"));
	topology = caprese.fit(data$genotypes,lambda,do.estimation);
	topology$data = data;
	#set rownames and colnames to the results
	rownames(topology$probabilities$marginal.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$marginal.probs) = "marginal probability";
	rownames(topology$probabilities$joint.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$joint.probs) = colnames(data$genotypes);
	rownames(topology$probabilities$conditional.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$conditional.probs) = "conditional probability";
	rownames(topology$parents.pos) = colnames(data$genotypes);
	colnames(topology$parents.pos) = "parent";
	rownames(topology$confidence) = c("probabilistic causality","hypergeometric test","adjusted hypergeometric test");
	colnames(topology$confidence) = "confidence";
	rownames(topology$confidence[[1,1]]) = colnames(data$genotypes);
	colnames(topology$confidence[[1,1]]) = colnames(data$genotypes);
	rownames(topology$confidence[[2,1]]) = colnames(data$genotypes);
	colnames(topology$confidence[[2,1]]) = colnames(data$genotypes);
	rownames(topology$confidence[[3,1]]) = colnames(data$genotypes);
	colnames(topology$confidence[[3,1]]) = colnames(data$genotypes);
	rownames(topology$adj.matrix) = colnames(data$genotypes);
	colnames(topology$adj.matrix) = colnames(data$genotypes);
	if(do.estimation==TRUE) {
		rownames(topology$probabilities$estimated.marginal.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$estimated.joint.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$estimated.joint.probs) = colnames(data$genotypes);
		rownames(topology$probabilities$estimated.conditional.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$estimated.conditional.probs) = "conditional probability";
	}
	#the reconstruction has been completed
	cat(paste("The reconstruction has been successfully completed.","\n"));
    return(topology);
}

#### end of file -- tronco.caprese.R
