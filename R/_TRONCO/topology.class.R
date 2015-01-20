#### topology.class.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @import methods
# The topology class definition
topology = NULL;
setClass("topology",
	representation(
		#the input valid dataset
		dataset = "data.frame",
		#the observed and estimated probabilities
		marginal.probs = "matrix",
		joint.probs = "matrix",
		conditional.probs = "matrix",
		estimated.marginal.probs = "matrix",
		estimated.joint.probs = "matrix",
		estimated.conditional.probs = "matrix",
		#the observed and estimated probabilities for the prima facie topology
		pf.marginal.probs = "matrix",
		pf.joint.probs = "matrix",
		pf.conditional.probs = "matrix",
		pf.estimated.marginal.probs = "matrix",
		pf.estimated.joint.probs = "matrix",
		pf.estimated.conditional.probs = "matrix",
		#the reconstructed topology
		adj.matrix = "matrix",
		parents.pos = "matrix",
		#the reconstructed prima facie topology
		pf.adj.matrix = "matrix",
		pf.parents.pos = "matrix",
		#the estimated error rates for the reconstructed and prima facie topologies
		error.fp = "numeric",
		error.fn = "numeric",
		pf.error.fp = "numeric",
		pf.error.fn = "numeric",
		#confidence in terms of prima facie scores
		confidence.scores = "list",
		#confidence in terms of non-parametric bootstrap
		confidence.np = "list",
		pf.confidence.np = "list",
		edge.confidence.np = "matrix",
		pf.edge.confidence.np = "matrix",
		bootstrap.settings.np = "list",
		bootstrap.np = "logical",
		#confidence in terms of parametric bootstrap
		confidence.p = "list",
		pf.confidence.p = "list",
		edge.confidence.p = "matrix",
		pf.edge.confidence.p = "matrix",
		bootstrap.settings.p = "list",
		bootstrap.p = "logical",
		#parameters
		parameters = "list",
		algorithm = "character"))

# A summary of the topology is displayed for this object
setMethod("show", "topology",
	function(object) {
		topology <- c();
		adj.matrix <- object@adj.matrix;
		names <- colnames(adj.matrix);
		for(i in 1:nrow(adj.matrix)) {
			for(j in 1:ncol(adj.matrix)) {
				if(adj.matrix[i,j] == 1) {
					topology <- c(topology,paste(names[i]," -> ",names[j],"\n"));
				}
			}
		}
		cat(object@algorithm,"progression model of",ncol(adj.matrix),"events.\n");
		cat("\n ");
		cat(topology);
		if(object@algorithm == "CAPRESE") {
			cat("\n Estimated false positives error rate:", object@error.fp);
			cat("\n Estimated false negative error rate:", object@error.fn);
		}
        else if(object@algorithm == "CAPRI") {
			cat("\n Estimated false positives \"bic\"error rate:", object@error.fp);
			cat("\n Estimated false negative \"bic\" error rate:", object@error.fp);
			cat("\n Estimated false positives \"prima facie\" error rate:", object@pf.error.fp);
			cat("\n Estimated false negative \"prima facie\" error rate:", object@pf.error.fn);
        }
	})

#### end of file -- topology.class.R
