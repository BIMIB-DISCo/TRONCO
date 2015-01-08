#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#verify the dataset and reconstruct the best dag topology
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#do.boot: should I perform bootstrap? Yes if TRUE, no otherwise
#nboot: integer number (greater than 0) of bootstrap sampling to be performed
#pvalue: pvalue for the tests (value between 0 and 1)
#do.estimation: should I perform the estimation of the error rates and probabilities?
#verbose: should I print the warnings? Yes if TRUE, no otherwise
#RETURN:
#topology: the reconstructed tree topology
"capri.fit" <-
function(dataset,do.boot,nboot,pvalue,do.estimation,verbose) {
	#reconstruct the causal topology
    #should I perform bootstrap? Yes if TRUE, no otherwise
    if(do.boot==TRUE) {
        prima.facie.parents = get.prima.facie.parents(dataset,nboot,pvalue);
    }
    else {
        prima.facie.parents = get.prima.facie.parents.no.boot(dataset);
    }
    #set the invalid entries in the adj.matrix
    if(exists("settings") && length(settings$data.values$num.events)>0) {
		if(nrow(prima.facie.parents$adj.matrix)-settings$data.values$num.events>0) {
			for (i in 1:nrow(prima.facie.parents$adj.matrix)) {
				for (j in (settings$data.values$num.events+1):nrow(prima.facie.parents$adj.matrix)) {
					prima.facie.parents$adj.matrix[i,j] = 0;
					if(do.boot==TRUE) {
						tmp = prima.facie.parents$pf.confidence[[1,1]];
						tmp[i,j] = 1;
						tmp[j,i] = 1;
						prima.facie.parents$pf.confidence[1,1] = list(tmp);
						tmp = prima.facie.parents$pf.confidence[[2,1]];
               			tmp[i,j] = 1;
               			tmp[j,i] = 1;
               			prima.facie.parents$pf.confidence[2,1] = list(tmp);
					}
				}
			}
		}
    }
    if(exists("settings") && length(settings$hypotheses$hlist)>0) {
		hypotheses.matrix = array(0,c(settings$hypotheses$num.hypotheses,settings$data.values$num.events));
		for (i in seq(1,(length(settings$hypotheses$hlist)-1),by=2)) {
			hypotheses.matrix[settings$hypotheses$hlist[i],settings$hypotheses$hlist[i+1]] = 1;
		}
		swap.matrix = prima.facie.parents$adj.matrix[(settings$data.values$num.events+1):nrow(prima.facie.parents$adj.matrix),1:settings$data.values$num.events];
		prima.facie.parents$adj.matrix[(settings$data.values$num.events+1):nrow(prima.facie.parents$adj.matrix),1:settings$data.values$num.events] = swap.matrix*hypotheses.matrix;
		for (i in 1:nrow(hypotheses.matrix)) {
			for(j in 1:nrow(hypotheses.matrix)) {
				if(hypotheses.matrix[i,j]==0) {
					if(do.boot==TRUE) {
               			tmp = prima.facie.parents$pf.confidence[[1,1]];
               			tmp[i+settings$data.values$num.events,j] = 1;
               			prima.facie.parents$pf.confidence[1,1] = list(tmp);
               			tmp = prima.facie.parents$pf.confidence[[2,1]];
               			tmp[i+settings$data.values$num.events,j] = 1;
               			prima.facie.parents$pf.confidence[2,1] = list(tmp);
					}
				}
			}
		}
	}
	#perform the likelihood fit by BIC score
	best.parents = perform.likelihood.fit(dataset,prima.facie.parents$adj.matrix);
	#set the conditional probabilities
	parents.pos.pf = array(list(),c(ncol(dataset),1));
    conditional.probs.pf = array(list(),c(ncol(dataset),1));
    parents.pos.bic = array(list(),c(ncol(dataset),1));
    conditional.probs.bic = array(list(),c(ncol(dataset),1));
    for(i in 1:ncol(dataset)) {
		for(j in 1:ncol(dataset)) {
			if(i!=j && best.parents$adj.matrix.prima.facie[i,j]==1) {
				parents.pos.pf[j,1] = list(c(unlist(parents.pos.pf[j,1]),i));
				conditional.probs.pf[j,1] = list(c(unlist(conditional.probs.pf[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
			if(i!=j && best.parents$adj.matrix.bic[i,j]==1) {
				parents.pos.bic[j,1] = list(c(unlist(parents.pos.bic[j,1]),i));
				conditional.probs.bic[j,1] = list(c(unlist(conditional.probs.bic[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
		}
    }
    parents.pos.pf[unlist(lapply(parents.pos.pf,is.null))] = list(-1);
    conditional.probs.pf[unlist(lapply(conditional.probs.pf,is.null))] = list(1);
    parents.pos.bic[unlist(lapply(parents.pos.bic,is.null))] = list(-1);
    conditional.probs.bic[unlist(lapply(conditional.probs.bic,is.null))] = list(1);
    if(do.estimation) {
		#estimate the error rates and, given them, the probabilities for the prima facie topology
		estimated.error.rates.pf = estimate.dag.error.rates(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.pf);
		estimated.probabilities.pf = estimate.dag.probs(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.pf,estimated.error.rates.pf);
		#estimate the error rates and, given them, the probabilities for the causal topology
		estimated.error.rates.bic = estimate.dag.error.rates(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.bic);
		estimated.probabilities.bic = estimate.dag.probs(dataset,prima.facie.parents$marginal.probs,prima.facie.parents$joint.probs,parents.pos.bic,estimated.error.rates.bic);
    }
    else {
		estimated.error.rates.pf = list(error.fp=NA,error.fn=NA);
		estimated.probabilities.pf = list(marginal.probs=NA,joint.probs=NA,conditional.probs=NA);
		estimated.error.rates.bic = list(error.fp=NA,error.fn=NA);
		estimated.probabilities.bic = list(marginal.probs=NA,joint.probs=NA,conditional.probs=NA);
    }
	#structures where to save the probabilities
    probabilities.pf = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.pf,estimated.marginal.probs=estimated.probabilities.pf$marginal.probs,estimated.joint.probs=estimated.probabilities.pf$joint.probs,estimated.conditional.probs=estimated.probabilities.pf$conditional.probs);
    probabilities.bic = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.bic,estimated.marginal.probs=estimated.probabilities.bic$marginal.probs,estimated.joint.probs=estimated.probabilities.bic$joint.probs,estimated.conditional.probs=estimated.probabilities.bic$conditional.probs);
    probabilities = list(probabilities.pf=probabilities.pf,probabilities.bic=probabilities.bic);
    parents.pos = list(parents.pos.pf=parents.pos.pf,parents.pos.bic=parents.pos.bic);
    error.rates = list(error.rates.pf=estimated.error.rates.pf,error.rates.bic=estimated.error.rates.bic);
	parameters = list(algorithm="CAPRI",do.boot=do.boot,nboot=nboot,pvalue=pvalue);
    #return the results
    topology = list(dataset=dataset,probabilities=probabilities,parents.pos=parents.pos,error.rates=error.rates,confidence=prima.facie.parents$pf.confidence,adj.matrix=best.parents,parameters=parameters);
	return(topology);
}

#### end of file -- capri.fit.R
