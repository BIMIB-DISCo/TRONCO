#### tronco.capri.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


tronco.capri <- function(dataset = NA, do.boot = TRUE, nboot = 100, pvalue = 0.05, do.estimation = TRUE, verbose = FALSE) {
	if(is.na(dataset)) {
		dataset = settings$data.values$dataset;
	}
	if(pvalue < 0 || pvalue > 1) {
		stop("The value of the pvalue has to be in [0:1]!",call.=FALSE);
	}
	#reconstruct the topology with CAPRI
	topology = capri.fit(dataset,do.boot,nboot,pvalue,do.estimation,verbose);
	#set the labels to the reconstructions
	#labels for the dataset
    colnames(topology$dataset) <- settings$visualization$labels;
    #labels for the reconstruction
    rownames(topology$adj.matrix$adj.matrix.prima.facie) <- settings$visualization$labels;
    colnames(topology$adj.matrix$adj.matrix.prima.facie) <- settings$visualization$labels;
    rownames(topology$parents.pos$parents.pos.pf) <- settings$visualization$labels;
    rownames(topology$adj.matrix$adj.matrix.bic) <- settings$visualization$labels;
    colnames(topology$adj.matrix$adj.matrix.bic) <- settings$visualization$labels;
    rownames(topology$parents.pos$parents.pos.bic) <- settings$visualization$labels;
    if(do.boot==TRUE) {
	    rownames(topology$confidence[[1,1]]) <- settings$visualization$labels;
		colnames(topology$confidence[[1,1]]) <- settings$visualization$labels;
	    rownames(topology$confidence[[2,1]]) <- settings$visualization$labels;
		colnames(topology$confidence[[2,1]]) <- settings$visualization$labels;
    }
    else {
	    	topology$confidence = array(list(array(NA, c(ncol(dataset),ncol(dataset)))), c(2,1));
    }
    #labels for the probabilities
    rownames(topology$probabilities$probabilities.pf$marginal.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$probabilities.pf$joint.probs) <- settings$visualization$labels;
    colnames(topology$probabilities$probabilities.pf$joint.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$probabilities.pf$conditional.probs) <- settings$visualization$labels;
    if(do.estimation==TRUE) {
		rownames(topology$probabilities$probabilities.pf$estimated.marginal.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$probabilities.pf$estimated.joint.probs) <- settings$visualization$labels;
		colnames(topology$probabilities$probabilities.pf$estimated.joint.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$probabilities.pf$estimated.conditional.probs) <- settings$visualization$labels;
    }
    else {
    		topology$probabilities$probabilities.pf$estimated.marginal.probs = matrix();
		topology$probabilities$probabilities.pf$estimated.joint.probs = matrix();
		topology$probabilities$probabilities.pf$estimated.conditional.probs = matrix();
		topology$error.rates$error.rates.bic = list(error.fp=numeric(),error.fn=numeric());
    }
    rownames(topology$probabilities$probabilities.bic$marginal.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$probabilities.bic$joint.probs) <- settings$visualization$labels;
    colnames(topology$probabilities$probabilities.bic$joint.probs) <- settings$visualization$labels;
    rownames(topology$probabilities$probabilities.bic$conditional.probs) <- settings$visualization$labels;
    if(do.estimation==TRUE) {
		rownames(topology$probabilities$probabilities.bic$estimated.marginal.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$probabilities.bic$estimated.joint.probs) <- settings$visualization$labels;
		colnames(topology$probabilities$probabilities.bic$estimated.joint.probs) <- settings$visualization$labels;
		rownames(topology$probabilities$probabilities.bic$estimated.conditional.probs) <- settings$visualization$labels;
    }
    else {
		topology$probabilities$probabilities.bic$estimated.marginal.probs = matrix();
		topology$probabilities$probabilities.bic$estimated.joint.probs = matrix();
		topology$probabilities$probabilities.bic$estimated.conditional.probs = matrix();
		topology$error.rates$error.rates.pf = list(error.fp=numeric(),error.fn=numeric());
    }
    #create the new object topology
    topology.obj <- new("topology",
    	#the input valid dataset
		dataset = topology$dataset,
		#the observed and estimated probabilities
		marginal.probs = topology$probabilities$probabilities.bic$marginal.probs,
		joint.probs = topology$probabilities$probabilities.bic$joint.probs,
		conditional.probs = topology$probabilities$probabilities.bic$conditional.probs,
		estimated.marginal.probs = topology$probabilities$probabilities.bic$estimated.marginal.probs,
		estimated.joint.probs = topology$probabilities$probabilities.bic$estimated.joint.probs,
		estimated.conditional.probs = topology$probabilities$probabilities.bic$estimated.conditional.probs,
		#the observed and estimated probabilities for the prima facie topology
		pf.marginal.probs = topology$probabilities$probabilities.pf$marginal.probs,
		pf.joint.probs = topology$probabilities$probabilities.pf$joint.probs,
		pf.conditional.probs = topology$probabilities$probabilities.pf$conditional.probs,
		pf.estimated.marginal.probs = topology$probabilities$probabilities.pf$estimated.marginal.probs,
		pf.estimated.joint.probs = topology$probabilities$probabilities.pf$estimated.joint.probs,
		pf.estimated.conditional.probs = topology$probabilities$probabilities.pf$estimated.conditional.probs,
		#the reconstructed topology
		adj.matrix = topology$adj.matrix$adj.matrix.bic,
		parents.pos = topology$parents.pos$parents.pos.bic,
		#the reconstructed prima facie topology
		pf.adj.matrix = topology$adj.matrix$adj.matrix.prima.facie,
		pf.parents.pos = topology$parents.pos$parents.pos.pf,
		#the estimated error rates for the reconstructed and prima facie topologies
		error.fp = topology$error.rates$error.rates.bic$error.fp,
		error.fn = topology$error.rates$error.rates.bic$error.fn,
		pf.error.fp = topology$error.rates$error.rates.pf$error.fp,
		pf.error.fn = topology$error.rates$error.rates.pf$error.fn,
		#confidence in terms of prima facie scores
		confidence.scores = list(temporal.priority=topology$confidence[[1,1]],probability.raising=topology$confidence[[2,1]]),
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
		algorithm = "CAPRI")
	#the reconstruction has been completed
	cat(paste("Executed CAPRI algorithm.\n"));
	if(do.estimation==TRUE) {
		cat("\nEstimated false positives \"bic\"error rate:",topology.obj@error.fp);
		cat("\nEstimated false negative \"bic\" error rate:",topology.obj@error.fn);
		cat("\nEstimated false positives \"prima facie\" error rate:",topology.obj@pf.error.fp);
		cat("\nEstimated false negative \"prima facie\" error rate:",topology.obj@pf.error.fn);
	}
    #save the reconstruction to the global workspace and return it
    assign("topology",topology.obj,envir=.GlobalEnv);
    return(topology.obj);
}

#### end of file -- tronco.capri.R
