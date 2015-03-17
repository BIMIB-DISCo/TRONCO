#### tronco.bootstrap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


tronco.bootstrap <- function( topology, type="non-parametric", nboot=100) {
	#check for the inputs to be correct
	if(is.null(topology)) {
		stop("A valid reconstruction has to be provided in order to estimate its confidence.",call.=FALSE);
    }
    if(topology$parameters$do.estimation==FALSE && type=="parametric") {
    		stop("To perform parametric bootstrap, the estimation of error rates and probabilities should be computed.",call.=FALSE);
    }
    #set all the needed parameters to perform the bootstrap
    if(type=="non-parametric" || type=="parametric") {
    		dataset = topology$data$genotypes;
		if(topology$parameters$algorithm=="CAPRESE") {
			lambda = topology$parameters$lambda;
			adj.matrix = topology$adj.matrix;
			if(type=="parametric") {
				if(topology$parameters$do.estimation==TRUE) {
					estimated.marginal.probs = topology$probabilities$estimated.marginal.probs;
					estimated.conditional.probs = topology$probabilities$estimated.conditional.probs;
					parents.pos = topology$parents.pos;
					error.rates = topology$error.rates;
				}
				else {
					stop("To perform parametric bootstrap, the estimation of the error rates should be performed first.",call.=FALSE);
				}
			}
		}
		else if(topology$parameters$algorithm=="CAPRI") {
			hypotheses = topology$data$hypotheses;
			if(is.null(hypotheses)) {
				hypotheses = NA;
			}
			command.capri = do.boot = topology$parameters$command;
	
			do.boot = topology$parameters$do.boot;
			nboot.capri = topology$parameters$nboot;
			pvalue = topology$parameters$pvalue;
			adj.matrix.pf = topology$adj.matrix$adj.matrix.pf;
			adj.matrix.bic = topology$adj.matrix$adj.matrix.bic;
			if(type=="parametric") {
				if(topology$parameters$do.estimation==TRUE) {
					estimated.marginal.probs.pf = topology$probabilities$probabilities.pf$estimated.marginal.probs;
					estimated.conditional.probs.pf = topology$probabilities$probabilities.pf$estimated.conditional.probs;
					parents.pos.pf = topology$parents.pos$parents.pos.pf;
					error.rates.pf = topology$error.rates$error.rates.pf;
					estimated.marginal.probs.bic = topology$probabilities$probabilities.bic$estimated.marginal.probs;
					estimated.conditional.probs.bic = topology$probabilities$probabilities.bic$estimated.conditional.probs;
					parents.pos.bic = topology$parents.pos$parents.pos.bic;
					error.rates.bic = topology$error.rates$error.rates.bic;
				}
				else {
					stop("To perform parametric bootstrap, the estimation of the error rates should be performed first.",call.=FALSE);
				}
			}
		}
    }
    else {
		stop("The types of bootstrap that can be performed are: non-parametric or parametric.",call.=FALSE);
    }
    #perform the selected bootstrap procedure
    cat("Executing now the bootstrap procedure, this may take a long time...\n");
    if(topology$parameters$algorithm=="CAPRESE") {
		if(type=="non-parametric") {
			curr.boot = bootstrap.caprese(dataset,lambda,adj.matrix,type,NA,NA,NA,nboot);
		}
		else if(type=="parametric") {
			curr.boot = bootstrap.caprese(dataset,lambda,adj.matrix,type,estimated.marginal.probs,estimated.conditional.probs,error.rates,nboot);
		}
		topology$bootstrap = curr.boot;
		cat("\nConfidence overall value: ",curr.boot$confidence$overall.value);
		cat("\nConfidence overall frequency: ",curr.boot$confidence$overall.frequency);
		cat(paste("\nPerformed ", type, " bootstrap with ",nboot," resampling and ",lambda," as shrinkage parameter.\n\n",sep =""));
    }
    else if(topology$parameters$algorithm=="CAPRI") {
		if(type=="non-parametric") {
			curr.boot = bootstrap.capri(dataset,hypotheses,command.capri,do.boot,nboot.capri,pvalue,adj.matrix.pf,adj.matrix.bic,type,NA,NA,NA,NA,NA,NA,NA,NA,nboot, 				REGULARIZATION=topology$parameters$REGULARIZATION);
		}
		else if(type=="parametric") {
			curr.boot = bootstrap.capri(dataset,hypotheses,command.capri,do.boot,nboot.capri,pvalue,adj.matrix.pf,adj.matrix.bic,type,estimated.marginal.probs.pf,estimated.conditional.probs.pf,parents.pos.pf,error.rates.pf,estimated.marginal.probs.bic,estimated.conditional.probs.bic,parents.pos.bic,error.rates.bic,nboot, REGULARIZATION=topology$parameters$REGULARIZATION);
		}
		topology$bootstrap = curr.boot;
		cat("\nConfidence overall \"prima facie\" value:",curr.boot$confidence$confidence.pf$overall.value.pf);
		cat("\nConfidence overall \"prima facie\" frequency:",curr.boot$confidence$confidence.pf$overall.frequency.pf);
		cat("\nConfidence overall \"bic\" value:",curr.boot$confidence$confidence.bic$overall.value.bic);
		cat("\nConfidence overall \"bic\" frequency:",curr.boot$confidence$confidence.bic$overall.frequency.bic);
		if(do.boot==TRUE) {
			cat(paste("\n\nPerformed ",type," bootstrap with ",nboot," resampling and ",pvalue," as pvalue for the statistical tests.\n\n",sep =""));
		}
		else {
			cat(paste("\n\nPerformed ",type," bootstrap with ",nboot," resampling.\n\n",sep =""));
		}
    }
    return(topology);
}

#### end of file -- tronco.bootstrap.R
