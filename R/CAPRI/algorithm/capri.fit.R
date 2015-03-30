#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# reconstruct the best dag topology running CAPRI algorithm
# INPUT:
# dataset: a dataset describing a progressive phenomenon
# hypotheses: hypotheses to be considered in the reconstruction
# command: type of search for the likelihood fit, either hill climbing (hc) or tabu (tabu)
# do.boot: should I perform bootstrap? Yes if TRUE, no otherwise
# nboot: integer number (greater than 0) of bootstrap sampling to be performed
# pvalue: pvalue for the tests (value between 0 and 1)
# min.boot: minimum number of bootstrapping to be performed
# min.stat: should I keep bootstrapping untill I have nboot valid values?
# boot.seed: seed to be used for the sampling
# do.estimation: should I perform the estimation of the error rates and probabilities?
# RETURN:
# topology: the reconstructed tree topology
"capri.fit" <-
function( dataset, hypotheses = NA, command = "hc", regularization = "bic", do.boot = TRUE, nboot = 100, pvalue = 0.05, min.boot = 3, min.stat = TRUE, boot.seed = 12345, do.estimation = FALSE ) {
	
	# structure with the set of valid edges
	# I start from the complete graph, i.e., I have no prior and all the connections are possibly causal
	adj.matrix = array(1,c(ncol(dataset),ncol(dataset)));
	colnames(adj.matrix) = colnames(dataset);
	rownames(adj.matrix) = colnames(dataset);
	
	# the diagonal of the adjacency matrix should not be considered, i.e., no self cause is allowed
	diag(adj.matrix) = 0;
	
	# consider any hypothesis
	adj.matrix = hypothesis.adj.matrix(hypotheses,adj.matrix);
	
	# reconstruct the prima facie topology
    # should I perform bootstrap? Yes if TRUE, no otherwise
    if(do.boot==TRUE) {
        prima.facie.parents = get.prima.facie.parents.do.boot(dataset,hypotheses,nboot,pvalue,adj.matrix,min.boot,min.stat,boot.seed);
    }
    else {
        prima.facie.parents = get.prima.facie.parents.no.boot(dataset,hypotheses,adj.matrix);
    }
    
	# perform the likelihood fit by BIC score on the prima facie topology
	best.parents = perform.likelihood.fit(dataset,prima.facie.parents$adj.matrix,command,regularization=regularization);
	
	# set the structure to save the conditional probabilities of the reconstructed topology
	parents.pos.pf = array(list(),c(ncol(dataset),1));
    conditional.probs.pf = array(list(),c(ncol(dataset),1));
    parents.pos.fit = array(list(),c(ncol(dataset),1));
    conditional.probs.fit = array(list(),c(ncol(dataset),1));
    
    # compute the conditional probabilities
    for(i in 1:ncol(dataset)) {
		for(j in 1:ncol(dataset)) {
			if(i!=j && best.parents$adj.matrix$adj.matrix.pf[i,j]==1) {
				parents.pos.pf[j,1] = list(c(unlist(parents.pos.pf[j,1]),i));
				conditional.probs.pf[j,1] = list(c(unlist(conditional.probs.pf[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
			if(i!=j && best.parents$adj.matrix$adj.matrix.fit[i,j]==1) {
				parents.pos.fit[j,1] = list(c(unlist(parents.pos.fit[j,1]),i));
				conditional.probs.fit[j,1] = list(c(unlist(conditional.probs.fit[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
		}
    }
    parents.pos.pf[unlist(lapply(parents.pos.pf,is.null))] = list(-1);
    conditional.probs.pf[unlist(lapply(conditional.probs.pf,is.null))] = list(1);
    parents.pos.fit[unlist(lapply(parents.pos.fit,is.null))] = list(-1);
    conditional.probs.fit[unlist(lapply(conditional.probs.fit,is.null))] = list(1);
    
    # perform the estimation of the probabilities if requested
    if(do.estimation) {
		# estimate the error rates and, given them, the probabilities for the prima facie topology
		estimated.error.rates.pf = estimate.dag.error.rates(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.pf);
		estimated.probabilities.pf = estimate.dag.probs(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.pf,estimated.error.rates.pf);
		# estimate the error rates and, given them, the probabilities for the causal topology
		estimated.error.rates.fit = estimate.dag.error.rates(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.fit);
		estimated.probabilities.fit = estimate.dag.probs(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.fit,estimated.error.rates.fit);
    }
    else {
		estimated.error.rates.pf = list(error.fp=NA,error.fn=NA);
		estimated.probabilities.pf = list(marginal.probs=NA,joint.probs=NA,conditional.probs=NA);
		estimated.error.rates.fit = list(error.fp=NA,error.fn=NA);
		estimated.probabilities.fit = list(marginal.probs=NA,joint.probs=NA,conditional.probs=NA);
    }
    
	# set structures where to save the results
    probabilities.pf = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.pf,estimated.marginal.probs=estimated.probabilities.pf$marginal.probs,estimated.joint.probs=estimated.probabilities.pf$joint.probs,estimated.conditional.probs=estimated.probabilities.pf$conditional.probs);
    probabilities.fit = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.fit,estimated.marginal.probs=estimated.probabilities.fit$marginal.probs,estimated.joint.probs=estimated.probabilities.fit$joint.probs,estimated.conditional.probs=estimated.probabilities.fit$conditional.probs);
    probabilities = list(probabilities.pf=probabilities.pf,probabilities.fit=probabilities.fit);
    parents.pos = list(parents.pos.pf=parents.pos.pf,parents.pos.fit=parents.pos.fit);
    error.rates = list(error.rates.pf=estimated.error.rates.pf,error.rates.fit=estimated.error.rates.fit);
	parameters = list(algorithm="CAPRI",command=command,regularization=regularization,do.boot=do.boot,nboot=nboot,pvalue=pvalue,min.boot=min.boot,min.stat=min.stat,boot.seed=boot.seed,do.estimation=do.estimation);
	
    # return the results
    topology = list(dataset=dataset,probabilities=probabilities,parents.pos=parents.pos,error.rates=error.rates,confidence=prima.facie.parents$pf.confidence,adj.matrix=best.parents$adj.matrix,parameters=parameters);
	return(topology);
	
}

#### end of file -- capri.fit.R
