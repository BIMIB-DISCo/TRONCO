#### tronco.capri.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


tronco.capri <- function(data, do.boot = TRUE, nboot = 100, pvalue = 0.05, do.estimation = FALSE) {
	#check for the inputs to be correct
	if(is.null(data) || is.null(data$genotypes)) {
		stop("The dataset given as input is not valid.");
	}
	if(is.null(data$hypotheses)) {
		data$hypotheses = NA;
	}
	if(pvalue < 0 || pvalue > 1) {
		stop("The value of the pvalue has to be in [0:1]!",call.=FALSE);
	}
	#reconstruct the topology with CAPRI
	cat(paste("Running CAPRI algorithm.","\n"));
	topology = capri.fit(data$genotypes,data$hypotheses,do.boot,nboot,pvalue,do.estimation);
	topology$data = data;
	#set rownames and colnames to the results
	rownames(topology$probabilities$probabilities.pf$marginal.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.pf$marginal.probs) = "marginal probability";
	rownames(topology$probabilities$probabilities.pf$joint.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.pf$joint.probs) = colnames(data$genotypes);
	rownames(topology$probabilities$probabilities.pf$conditional.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.pf$conditional.probs) = "conditional probability";
	rownames(topology$probabilities$probabilities.bic$marginal.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.bic$marginal.probs) = "marginal probability";
	rownames(topology$probabilities$probabilities.bic$joint.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.bic$joint.probs) = colnames(data$genotypes);
	rownames(topology$probabilities$probabilities.bic$conditional.probs) = colnames(data$genotypes);
	colnames(topology$probabilities$probabilities.bic$conditional.probs) = "conditional probability";
	rownames(topology$parents.pos$parents.pos.pf) = colnames(data$genotypes);
	colnames(topology$parents.pos$parents.pos.pf) = "parent";
	rownames(topology$parents.pos$parents.pos.bic) = colnames(data$genotypes);
	colnames(topology$parents.pos$parents.pos.bic) = "parent";
	rownames(topology$cpt$cpt.pf) = colnames(data$genotypes);
	colnames(topology$cpt$cpt.pf) = "parent";
	rownames(topology$cpt$cpt.bic) = colnames(data$genotypes);
	colnames(topology$cpt$cpt.bic) = "parent";
	rownames(topology$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
	colnames(topology$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
	rownames(topology$adj.matrix$adj.matrix.bic) = colnames(data$genotypes);
	colnames(topology$adj.matrix$adj.matrix.bic) = colnames(data$genotypes);
	rownames(topology$confidence) = c("temporal priority","probability raising","hypergeometric test","adjusted hypergeometric test");
	colnames(topology$confidence) = "confidence";
	if(do.boot==TRUE) {
		rownames(topology$confidence[[1,1]]) = colnames(data$genotypes);
		colnames(topology$confidence[[1,1]]) = colnames(data$genotypes);
		rownames(topology$confidence[[2,1]]) = colnames(data$genotypes);
		colnames(topology$confidence[[2,1]]) = colnames(data$genotypes);
	}
	rownames(topology$confidence[[3,1]]) = colnames(data$genotypes);
	colnames(topology$confidence[[3,1]]) = colnames(data$genotypes);
	rownames(topology$confidence[[4,1]]) = colnames(data$genotypes);
	colnames(topology$confidence[[4,1]]) = colnames(data$genotypes);
	if(do.estimation==TRUE) {
		rownames(topology$probabilities$probabilities.pf$estimated.marginal.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$probabilities.pf$estimated.joint.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.joint.probs) = colnames(data$genotypes);
		rownames(topology$probabilities$probabilities.pf$estimated.conditional.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.conditional.probs) = "conditional probability";
		rownames(topology$probabilities$probabilities.bic$estimated.marginal.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$probabilities.bic$estimated.joint.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.joint.probs) = colnames(data$genotypes);
		rownames(topology$probabilities$probabilities.bic$estimated.conditional.probs) = colnames(data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.conditional.probs) = "conditional probability";
	}
	#the reconstruction has been completed
	cat(paste("The reconstruction has been successfully completed.","\n"));
    return(topology);
}

#### end of file -- tronco.capri.R
