#### tronco.estimation.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


tronco.estimation <- function( topology, error.rates = NA ) {
	#check for the inputs to be correct
	if(is.null(topology)) {
		stop("A valid reconstruction has to be provided in order to estimate its confidence.",call.=FALSE);
    }
    #run the estimations for the required algorithm
    if(topology$parameters$algorithm=="CAPRESE") {
    		#if I also need to estimate the error rates
    		if(is.na(error.rates[1])) {
    			#estimate the error rates
    			error.rates = estimate.tree.error.rates(topology$probabilities$marginal.probs,topology$probabilities$joint.probs,topology$parents.pos);
    		}
    		#estimate the probabilities given the error rates
		estimated.probabilities = estimate.tree.probs(topology$probabilities$marginal.probs,topology$probabilities$joint.probs,topology$parents.pos,error.rates);
		#set the estimated error rates and probabilities
		topology$error.rates = error.rates;
		topology$probabilities$estimated.marginal.probs = estimated.probabilities$marginal.probs;
		topology$probabilities$estimated.joint.probs = estimated.probabilities$joint.probs;
		topology$probabilities$estimated.conditional.probs = estimated.probabilities$conditional.probs;
		#set colnames and rownames
		rownames(topology$probabilities$estimated.marginal.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$estimated.joint.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$estimated.joint.probs) = colnames(topology$data$genotypes);
		rownames(topology$probabilities$estimated.conditional.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$estimated.conditional.probs) = "conditional probability";
    }
    else if(topology$parameters$algorithm=="CAPRI") {
    		#if I also need to estimate the error rates
    		if(is.na(error.rates[1])) {
    			#estimate the error rates
    			estimated.error.rates.pf = estimate.dag.error.rates(topology$data$genotypes,topology$probabilities$probabilities.pf$marginal.probs,topology$probabilities$probabilities.pf$joint.probs,topology$parents.pos$parents.pos.pf);
    			estimated.error.rates.bic = estimate.dag.error.rates(topology$data$genotypes,topology$probabilities$probabilities.bic$marginal.probs,topology$probabilities$probabilities.bic$joint.probs,topology$parents.pos$parents.pos.bic);
    			error.rates = list(error.rates.pf=estimated.error.rates.pf,error.rates.bic=estimated.error.rates.bic);
    		}
    		#estimate the probabilities given the error rates
		estimated.probabilities.pf = estimate.dag.probs(topology$data$genotypes,topology$probabilities$probabilities.pf$marginal.probs,topology$probabilities$probabilities.pf$joint.probs,topology$parents.pos$parents.pos.pf,error.rates$error.rates.pf);
		estimated.probabilities.bic = estimate.dag.probs(topology$data$genotypes,topology$probabilities$probabilities.bic$marginal.probs,topology$probabilities$probabilities.bic$joint.probs,topology$parents.pos$parents.pos.bic,error.rates$error.rates.bic);
		#set the estimated error rates and probabilities
		topology$error.rates = error.rates;
		topology$probabilities$probabilities.pf$estimated.marginal.probs = estimated.probabilities.pf$marginal.probs;
		topology$probabilities$probabilities.pf$estimated.joint.probs = estimated.probabilities.pf$joint.probs;
		topology$probabilities$probabilities.pf$estimated.conditional.probs = estimated.probabilities.bic$conditional.probs;
		topology$probabilities$probabilities.bic$estimated.marginal.probs = estimated.probabilities.bic$marginal.probs;
		topology$probabilities$probabilities.bic$estimated.joint.probs = estimated.probabilities.bic$joint.probs;
		topology$probabilities$probabilities.bic$estimated.conditional.probs = estimated.probabilities.bic$conditional.probs;
		#set colnames and rownames
		rownames(topology$probabilities$probabilities.pf$estimated.marginal.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$probabilities.pf$estimated.joint.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.joint.probs) = colnames(topology$data$genotypes);
		rownames(topology$probabilities$probabilities.pf$estimated.conditional.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.pf$estimated.conditional.probs) = "conditional probability";
		rownames(topology$probabilities$probabilities.bic$estimated.marginal.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.marginal.probs) = "marginal probability";
		rownames(topology$probabilities$probabilities.bic$estimated.joint.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.joint.probs) = colnames(topology$data$genotypes);
		rownames(topology$probabilities$probabilities.bic$estimated.conditional.probs) = colnames(topology$data$genotypes);
		colnames(topology$probabilities$probabilities.bic$estimated.conditional.probs) = "conditional probability";
    }
    else {
    		stop("A valid algorithm has to be provided in order to estimate its confidence.",call.=FALSE);
    }
    topology$parameters$do.estimation = TRUE;
    return(topology);
}

#### end of file -- tronco.estimation.R
