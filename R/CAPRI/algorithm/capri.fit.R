#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#reconstruct the best dag topology running CAPRI algorithm
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#hypotheses: hypotheses to be considered in the reconstruction
#command: type of search for the likelihood fit, either hill climbing (hc) or tabu (tabu)
#do.boot: should I perform bootstrap? Yes if TRUE, no otherwise
#nboot: integer number (greater than 0) of bootstrap sampling to be performed
#pvalue: pvalue for the tests (value between 0 and 1)
#do.estimation: should I perform the estimation of the error rates and probabilities?
#RETURN:
#topology: the reconstructed tree topology
"capri.fit" <-
function(dataset, hypotheses = NA, command = "hc", do.boot = TRUE, nboot = 100, pvalue = 0.05, do.estimation = FALSE) {
	
	#structure with the set of valid edges
	#I start from the complete graph, i.e., I have no prior and all the connections are possibly causal
	adj.matrix = array(1,c(ncol(dataset),ncol(dataset)));
	#the diagonal of the adjacency matrix should not be considered, i.e., no self cause is allowed
	diag(adj.matrix) = 0;
	#consider the hypotheses if any
	if(!is.na(hypotheses[1])) {
		#set the invalid entries in the adj.matrix
		#hypotheses can not be causing other hypotheses
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix),(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix)] = 0;
		#consider the given hypotheses only toward the specified possible effects
		hypotheses.matrix = array(0,c(hypotheses$num.hypotheses,ncol(adj.matrix)-hypotheses$num.hypotheses));		
		for (i in 1:nrow(hypotheses$hlist)) {
			cause = which(hypotheses$hlist[i,1]==colnames(dataset));
			effect = which(hypotheses$hlist[i,2]==colnames(dataset));
			if(length(cause)>0 && length(effect)>0) {
				hypotheses.matrix[cause-ncol(dataset)+hypotheses$num.hypotheses,effect] = 1;
			}
		}
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix),1:(ncol(adj.matrix)-hypotheses$num.hypotheses)] = hypotheses.matrix;
		for(j in (ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix)) {
			for(k in 1:(ncol(adj.matrix)-hypotheses$num.hypotheses)) {
				if(adj.matrix[j,k] == 0) {
					adj.matrix[k,j] = 0;
				}
			}
		}
	}
	
	#reconstruct the prima facie topology
    #should I perform bootstrap? Yes if TRUE, no otherwise
    if(do.boot==TRUE) {
        prima.facie.parents = get.prima.facie.parents.do.boot(dataset,nboot,pvalue,adj.matrix);
    }
    else {
        prima.facie.parents = get.prima.facie.parents.no.boot(dataset,adj.matrix);
    }
    
	#perform the likelihood fit by BIC score on the prima facie topology
	best.parents = perform.likelihood.fit(dataset,prima.facie.parents$adj.matrix,command);
	
	#set the structure to save the conditional probabilities and the confidences of the reconstructed topology
	parents.pos.pf = array(list(),c(ncol(dataset),1));
    conditional.probs.pf = array(list(),c(ncol(dataset),1));
    parents.pos.bic = array(list(),c(ncol(dataset),1));
    conditional.probs.bic = array(list(),c(ncol(dataset),1));
    confidence <- array(list(), c(4,1));
    if(do.boot==TRUE) {
    		confidence[[1,1]] = prima.facie.parents$pf.confidence[[1,1]];
    		confidence[[2,1]] = prima.facie.parents$pf.confidence[[2,1]];
    }
    else {
    		confidence[[1,1]] = NA;
    		confidence[[2,1]] = NA;
    }
    confidence[[3,1]] = array(0, dim=c(ncol(dataset),ncol(dataset)));
    confidence[[4,1]] = array(0, dim=c(ncol(dataset),ncol(dataset)));
    
    #compute the conditional probabilities and the confidences
    hypergeometric.pvalues = vector();
    for(i in 1:ncol(dataset)) {
		for(j in 1:ncol(dataset)) {
			if(i!=j && best.parents$adj.matrix$adj.matrix.pf[i,j]==1) {
				parents.pos.pf[j,1] = list(c(unlist(parents.pos.pf[j,1]),i));
				conditional.probs.pf[j,1] = list(c(unlist(conditional.probs.pf[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
			if(i!=j && best.parents$adj.matrix$adj.matrix.bic[i,j]==1) {
				parents.pos.bic[j,1] = list(c(unlist(parents.pos.bic[j,1]),i));
				conditional.probs.bic[j,1] = list(c(unlist(conditional.probs.bic[j,1]),prima.facie.parents$joint.probs[i,j]/prima.facie.parents$marginal.probs[i]));
			}
        		if(i<j) {
        			#compute the confidence by hypergeometric test
        			confidence[[3,1]][i,j] = phyper(prima.facie.parents$joint.probs[i,j]*nrow(dataset),prima.facie.parents$marginal.probs[i]*nrow(dataset),nrow(dataset)-prima.facie.parents$marginal.probs[i]*nrow(dataset),prima.facie.parents$marginal.probs[j]*nrow(dataset),lower.tail=FALSE);
        			confidence[[3,1]][j,i] = confidence[[3,1]][i,j];
        			#save all the valid pvalues
        			hypergeometric.pvalues = append(hypergeometric.pvalues,confidence[[3,1]][i,j]);
        		}
		}
    }    
    parents.pos.pf[unlist(lapply(parents.pos.pf,is.null))] = list(-1);
    conditional.probs.pf[unlist(lapply(conditional.probs.pf,is.null))] = list(1);
    parents.pos.bic[unlist(lapply(parents.pos.bic,is.null))] = list(-1);
    conditional.probs.bic[unlist(lapply(conditional.probs.bic,is.null))] = list(1);
    
    #perform false discovery rate on the valid pvalues to compute the adjusted confidence
    hypergeometric.pvalues = p.adjust(hypergeometric.pvalues,method="fdr");
    #save the resulting pvalues
    cont = 0;
    for(i in 1:ncol(dataset)) {
        for (j in i:ncol(dataset)) {
        		if(i!=j) {
        			cont = cont + 1;
        			confidence[[4,1]][i,j] = hypergeometric.pvalues[cont];
        			confidence[[4,1]][j,i] = confidence[[4,1]][i,j];
        		}
        }
    }
    
    #perform the estimation of the probabilities if requested
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
    
	#set structures where to save the probabilities and the other results
    probabilities.pf = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.pf,estimated.marginal.probs=estimated.probabilities.pf$marginal.probs,estimated.joint.probs=estimated.probabilities.pf$joint.probs,estimated.conditional.probs=estimated.probabilities.pf$conditional.probs);
    probabilities.bic = list(marginal.probs=prima.facie.parents$marginal.probs,joint.probs=prima.facie.parents$joint.probs,conditional.probs=conditional.probs.bic,estimated.marginal.probs=estimated.probabilities.bic$marginal.probs,estimated.joint.probs=estimated.probabilities.bic$joint.probs,estimated.conditional.probs=estimated.probabilities.bic$conditional.probs);
    probabilities = list(probabilities.pf=probabilities.pf,probabilities.bic=probabilities.bic);
    parents.pos = list(parents.pos.pf=parents.pos.pf,parents.pos.bic=parents.pos.bic);
    error.rates = list(error.rates.pf=estimated.error.rates.pf,error.rates.bic=estimated.error.rates.bic);
	parameters = list(algorithm="CAPRI",command=command,do.boot=do.boot,nboot=nboot,pvalue=pvalue,do.estimation=do.estimation);
	
    #return the results
    topology = list(data=dataset,probabilities=probabilities,parents.pos=parents.pos,error.rates=error.rates,confidence=confidence,adj.matrix=best.parents$adj.matrix,parameters=parameters);
	return(topology);
	
}

#### end of file -- capri.fit.R
