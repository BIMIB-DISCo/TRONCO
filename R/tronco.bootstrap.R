#### tronco.bootstrap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export tronco.bootstrap
#' @title perform bootstrap estimation
#'
#' @description
#' \code{tronco.bootstrap} perform parametric and non-parametric bootstrap procedures 
#' 
#' @param curr.topology A reconstructed curr.topology
#' @param type The type of bootstrap to be performed, non-parametric and parametric types are implemented.
#' @param nboot Number of bootstrap resampling.
#' @return A curr.topology object with bootstrap estimation
#' 
tronco.bootstrap <- function(curr.topology=NA,dataset=NA,lambda=NA,do.boot=NA,nboot.capri=NA,pvalue=NA,reconstructed.curr.topology.pf=NA,reconstructed.curr.topology=NA,type="non-parametric",estimated.marginal.probs.pf=NA,estimated.conditional.probs.pf=NA,parents.pos.pf=NA,error.rates.pf=NA,estimated.marginal.probs=NA, estimated.conditional.probs=NA,parents.pos=NA,error.rates=NA,nboot=100) {
	if(suppressWarnings(is.na(curr.topology))) {
		curr.topology = topology;
	}
	if(is.null(curr.topology)) {
		stop("Missing parameters for the function tronco.bootstrap",call.=FALSE);
    }
    #set all the needed parameters
    if(type=="non-parametric" || type=="parametric") {
		if(is.na(dataset)) {
			dataset = curr.topology@dataset;
		}
		if(is.na(reconstructed.curr.topology)) {
			reconstructed.curr.topology = curr.topology@adj.matrix;
		}
		if(curr.topology@algorithm=="CAPRESE") {
			if(is.na(lambda)) {
				lambda = curr.topology@parameters$lambda;
			}
			else {
				if(lambda < 0 || lambda > 1) {
					stop("The value of lambda has to be in [0:1]!",call.=FALSE);
				}
			}
			if(type=="parametric") {
				if(is.na(estimated.marginal.probs)) {
					estimated.marginal.probs = curr.topology@marginal.probs;
				}
				if(is.na(estimated.conditional.probs)) {
					estimated.conditional.probs = curr.topology@conditional.probs;
				}
				if(is.na(parents.pos)) {
					parents.pos = curr.topology@parents.pos;
				}
				if(is.na(error.rates)) {
					error.rates = list(error.fp=curr.topology@error.fp,error.fn=curr.topology@error.fn);
				}
			}
		}
		else if(curr.topology@algorithm=="CAPRI") {
			if(is.na(do.boot)) {
				do.boot = curr.topology@parameters$do.boot;
			}
			if(is.na(nboot.capri)) {
			nboot.capri = curr.topology@parameters$nboot;
			}
			if(is.na(pvalue)) {
				pvalue = curr.topology@parameters$pvalue;
			}
			else {
				if(pvalue < 0 || lambda > 1) {
					stop("The value of the pvalue has to be in [0:1]!",call.=FALSE);
				}
			}
			if(is.na(reconstructed.curr.topology.pf)) {
				reconstructed.curr.topology.pf = curr.topology@pf.adj.matrix;;
			}
			if(type=="parametric") {
				if(is.na(estimated.marginal.probs)) {
					estimated.marginal.probs = curr.topology@marginal.probs;
				}
				if(is.na(estimated.conditional.probs)) {
					estimated.conditional.probs = curr.topology@conditional.probs;
				}
				if(is.na(parents.pos)) {
					parents.pos = curr.topology@parents.pos;
				}
				if(is.na(error.rates)) {
					error.rates = list(error.fp=curr.topology@error.fp,error.fn=curr.topology@error.fn);
				}
				if(is.na(estimated.marginal.probs.pf)) {
					estimated.marginal.probs.pf = curr.topology@pf.marginal.probs;
				}
				if(is.na(estimated.conditional.probs.pf)) {
					estimated.conditional.probs.pf = curr.topology@pf.conditional.probs;
				}
				if(is.na(parents.pos.pf)) {
					parents.pos.pf = curr.topology@pf.parents.pos;
				}
				if(is.na(error.rates.pf)) {
					error.rates.pf = list(error.fp=curr.topology@pf.error.fp,error.fn=curr.topology@pf.error.fn);
				}
			}
		}
    }
    else {
		stop("The valid types of bootstrap are: non-parametric and parametric.",call.=FALSE);
    }
    #perform the selected bootstrap procedure
    cat("Executing now the bootstrap procedure, this may take a long time...\n");
    if(curr.topology@algorithm=="CAPRESE") {
		curr.boot = bootstrap.caprese(dataset,lambda,reconstructed.curr.topology,type,estimated.marginal.probs,estimated.conditional.probs,error.rates,nboot);
		if(type=="non-parametric") {
			curr.topology@confidence.np <- curr.boot$confidence;
			curr.topology@edge.confidence.np <- curr.boot$edge.confidence;
			curr.topology@bootstrap.settings.np <- curr.boot$bootstrap.settings;
			curr.topology@bootstrap.np <-  TRUE;
		}
		else if(type=="parametric") {
			curr.topology@confidence.p <- curr.boot$confidence;
			curr.topology@edge.confidence.p <- curr.boot$edge.confidence;
			curr.topology@bootstrap.settings.p <- curr.boot$bootstrap.settings;
			curr.topology@bootstrap.p <-  TRUE;
		}
		cat("\nConfidence overall value:",curr.boot$confidence$overall.value);
		cat("\nConfidence overall frequency:",curr.boot$confidence$overall.frequency);
		cat(paste("\n\nExecuted ", type, " bootstrap with ",nboot," resampling and ",lambda," as shrinkage parameter.\n\n",sep =""));
    }
    else if(curr.topology@algorithm=="CAPRI") {
		curr.boot = bootstrap.capri(dataset,do.boot,nboot.capri,pvalue,reconstructed.curr.topology.pf,reconstructed.curr.topology,type,estimated.marginal.probs.pf,estimated.conditional.probs.pf,parents.pos.pf,error.rates.pf,estimated.marginal.probs,estimated.conditional.probs,parents.pos,error.rates,nboot);
		if(type=="non-parametric") {
			curr.topology@confidence.np <- curr.boot$confidence$confidence.bic;
			curr.topology@pf.confidence.np <- curr.boot$confidence$confidence.pf;
			curr.topology@edge.confidence.np <- curr.boot$edge.confidence$edge.confidence.bic;
			curr.topology@pf.edge.confidence.np <- curr.boot$edge.confidence$edge.confidence.pf;
			curr.topology@bootstrap.settings.np <- curr.boot$bootstrap.settings;
			curr.topology@bootstrap.np <- TRUE;
		}
		else if(type=="parametric") {
			curr.topology@confidence.p <- curr.boot$confidence$confidence.bic;
			curr.topology@pf.confidence.p <- curr.boot$confidence$confidence.pf;
			curr.topology@edge.confidence.p <- curr.boot$edge.confidence$edge.confidence.bic;
			curr.topology@pf.edge.confidence.p <- curr.boot$edge.confidence$edge.confidence.pf;
			curr.topology@bootstrap.settings.p <- curr.boot$bootstrap.settings;
			curr.topology@bootstrap.p <- TRUE;
		}
		cat("\nConfidence overall \"bic\" value:",curr.boot$confidence$confidence.bic$overall.value.bic);
		cat("\nConfidence overall \"bic\" frequency:",curr.boot$confidence$confidence.bic$overall.frequency.bic);
		cat("\nConfidence overall \"prima facie\" value:",curr.boot$confidence$confidence.pf$overall.value.pf);
		cat("\nConfidence overall \"prima facie\" frequency:",curr.boot$confidence$confidence.pf$overall.frequency.pf);
		if(do.boot==TRUE) {
			cat(paste("\n\nExecuted ",type," bootstrap with ",nboot," resampling and ",pvalue," as pvalue for the statistical tests.\n\n",sep =""));
		}
		else {
			cat(paste("\n\nExecuted ",type," bootstrap with ",nboot," resampling.\n\n",sep =""));
		}
    }
    #save the reconstruction to the global workspace and return it
    assign("curr.topology", curr.topology,envir=.GlobalEnv);
    return(curr.topology);
}

#### end of file -- tronco.bootstrap.R
