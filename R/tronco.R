#' Genotype-level cancer progression models describe the ordering of
#' accumulating mutations, e.g., somatic mutations / copy number variations,
#' during cancer development. These graphical models help understand the
#' causal structure involving events promoting cancer progression, possibly
#' predicting complex patterns characterising genomic progression of a cancer.
#' Reconstructed models can be used to better characterise genotype-phenotype
#' relation, and suggest novel targets for therapy design. 
#'
#' TRONCO (TRanslational ONCOlogy) is a R package aimed at collecting
#' state-of-the-art algorithms to infer progression models from
#' cross-sectional data, i.e., data collected from independent patients which
#' does not necessarily incorporate any evident temporal information. These
#' algorithms require a binary input matrix where: (i) each row represents a
#' patient genome, (ii) each column an event relevant to the progression (a
#' priori selected) and a 0/1 value models the absence/presence of a certain
#' mutation in a certain patient. The current first version of TRONCO
#' implements the CAPRESE algorithm (Cancer PRogression Extraction with Single
#' Edges) to infer possible progression models arranged as trees; cfr.
#' Inferring tree causal models of cancer progression with probability
#' raising, L. Olde Loohuis, G. Caravagna, A. Graudenzi, D. Ramazzotti, G.
#' Mauri, M. Antoniotti and B. Mishra. PLoS One, to appear. This vignette
#' shows how to use TRONCO to infer a tree model of ovarian cancer progression
#' from CGH data of copy number alterations (classified as gains or losses
#' over chromosome's arms). The dataset used is available in the SKY/M-FISH
#' database.
#'
#' @docType package
#' @name TRONCO
NULL

#' Reconstruc a progression model using CAPRESE algorithm
#'
#' @examples
#' data(maf)
#' mutations = import.MAF(maf)
#' recon = tronco.caprese(mutations)
#' tronco.plot(recon)
#'
#' @title tronco caprese
#' @param data A TRONCO compliant dataset.
#' @param lambda TODO
#' @param do.estimation TODO
#' @param silent TODO
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.caprese
tronco.caprese <- function( data, lambda = 0.5, do.estimation = FALSE, silent = FALSE ) {

	#check for the inputs to be correct
    if(is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    if(lambda < 0 || lambda > 1) {
        stop("The value of the shrinkage parameter lambda has to be in [0:1]!",call.=FALSE);
    }
    
    # check for the input to be compliant
    is.compliant(data)
    
    #reconstruct the reconstruction with CAPRESE
    if(silent==FALSE) {
	    cat(paste0(
	        '*** Inferring a progression model with the following settings.\n',
	        '\tDataset size: n = ', nsamples(data), ', m = ', nevents(data), '.\n',
	        '\tAlgorithm: CAPRESE with shrinkage coefficient: ', lambda, '.\n'
	    ))
	}
    reconstruction = caprese.fit(data$genotypes,lambda,do.estimation,silent);
    
    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    
    for (i in 1:length(reconstruction$model)) {
        
        #set rownames and colnames to the probabilities
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";
        
        #set rownames and colnames to the parents positions
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";
        
        #set rownames and colnames to the adjacency matrices
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        
        if(do.estimation==TRUE) {
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.marginal.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.marginal.probs) = "marginal probability";
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.conditional.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.conditional.probs) = "conditional probability";
        }
        
    }
    
    # structure to save the results
    results = data;
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;
    
    # the reconstruction has been completed
	if(!silent) cat(paste(
	  "The reconstruction has been successfully completed in ", 
	  round(reconstruction$execution.time[3]/60,digits=0), 
	  "minutes.\n"));
	
  return(results);
    
}



#' Reconstruc a progression model using CAPRI algorithm
#'
#' @examples
#' data(maf)
#' mutations = import.MAF(maf)
#' recon = tronco.capri(mutations)
#' tronco.plot(recon)
#'
#' @title tronco capri
#' @param data A TRONCO compliant dataset.
#' @param command TODO
#' @param regularization TODO
#' @param do.boot TODO
#' @param nboot TODO
#' @param pvalue TODO
#' @param min.boot TODO
#' @param min.stat TODO
#' @param boot.seed TODO
#' @param do.estimation TODO
#' @param silent TODO
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.capri
tronco.capri <- function( data, 
    command = "hc", 
    regularization = c("bic","aic"), 
    do.boot = TRUE, 
    nboot = 100, 
    pvalue = 0.05, 
    min.boot = 3, 
    min.stat = TRUE, 
    boot.seed = NULL, 
    do.estimation = FALSE, 
    silent = FALSE ) 
{
    #check for the inputs to be correct
    if(is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    if(is.null(data$hypotheses)) {
        data$hypotheses = NA;
    }
    if(command != "hc" && command != "tabu") {
        stop("The inference can be performed either by hill climbing or tabu search!",call.=FALSE);
    }
    if(pvalue < 0 || pvalue > 1) {
        stop("The value of the pvalue has to be in [0:1]!",call.=FALSE);
    }
    
    # check for the input to be compliant
    is.compliant(data)
    
    # reconstruct the reconstruction with CAPRI
    if(is.null(boot.seed)) {
        my.seed = "NULL"    
    }
    else {
            my.seed = boot.seed;
    }
    if(silent==FALSE) {
	    cat(paste0(
	        '*** Inferring a progression model with the following settings.\n',
	        '\tDataset size: n = ', nsamples(data), ', m = ', nevents(data), '.\n',
	        '\tAlgorithm: CAPRI with \"', paste0(regularization,collapse=", "), '\" regularization and \"', command, '\" likelihood-fit strategy.\n',
	        '\tRandom seed: ', my.seed, '.\n',
	        '\tBootstrap iterations (Wilcoxon): ', ifelse(do.boot, nboot, 'disabled'), '.\n',
	        ifelse(do.boot, 
	            paste0('\t\texhaustive bootstrap: ', min.stat, '.\n\t\tp-value: ', pvalue, '.\n\t\tminimum bootstrapped scores: ', min.boot, '.\n'), '')        
	        ))
	}
        
    reconstruction = capri.fit(data$genotypes,data$hypotheses,command=command,regularization=regularization,do.boot=do.boot,nboot=nboot,pvalue=pvalue,min.boot=min.boot,min.stat=min.stat,boot.seed=boot.seed,do.estimation=do.estimation,silent=silent);
    
    # reconstruction[[1]] = data;
    # names(reconstruction)[1] = "data"
    # reconstruction$hypotheses = NULL;
    
    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    
    for (i in 1:length(reconstruction$model)) {
        
        #set rownames and colnames to the probabilities
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";
        
        #set rownames and colnames to the parents positions
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";
        
        #set rownames and colnames to the adjacency matrices
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        
        if(do.estimation==TRUE) {
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.marginal.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.marginal.probs) = "marginal probability";
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            rownames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.conditional.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[i]]$probabilities$probabilities.fit$estimated.conditional.probs) = "conditional probability";
        }
        
    }
    
    # structure to save the results
    results = data;
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;
    
  

	
  
    # the reconstruction has been completed
    if(!silent) cat(paste(
      "The reconstruction has been successfully completed in ", 
      round(reconstruction$execution.time[3]/60,digits=0), 
      "minutes.\n"));
  
    return(results);
}

#### end of file -- tronco.capri.R


#### tronco.estimation.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

#' @export
tronco.estimation <- function( reconstruction, error.rates = NA ) {
	
    # check for the inputs to be correct
    if(is.null(reconstruction)) {
        stop("A valid reconstruction has to be provided in order to estimate its confidence.",call.=FALSE);
    }
    
    # check for the input to be compliant
    is.compliant(reconstruction)
    
    #run the estimations for the required algorithm
    if(reconstruction$parameters$algorithm=="CAPRESE") {
    	
    		cat("Executing now the estimation procedure, this may take a long time...\n")
    	
    		# if I also need to estimate the error rates
    		if(is.na(error.rates[1])) {
    			# estimate the error rates
    			error.rates = estimate.tree.error.rates(as.marginal.probs(reconstruction,models="caprese")[[1]],as.joint.probs(reconstruction,models="caprese")[[1]],as.parents.pos(reconstruction,models="caprese")[[1]]);
    		}
    		
    		# estimate the probabilities given the error rates
        estimated.probabilities = estimate.tree.probs(as.marginal.probs(reconstruction,models="caprese")[[1]],as.joint.probs(reconstruction,models="caprese")[[1]],as.parents.pos(reconstruction,models="caprese")[[1]],error.rates);
        
        # set the estimated error rates and probabilities
        probabilities.fit = list(estimated.marginal.probs=estimated.probabilities$marginal.probs,estimated.joint.probs=estimated.probabilities$joint.probs,estimated.conditional.probs=estimated.probabilities$conditional.probs);
        
        reconstruction$model[["caprese"]]$error.rates = error.rates
        reconstruction$model[["caprese"]]$probabilities$probabilities.fit = probabilities.fit
        
        # set colnames and rownames
       	rownames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.marginal.probs) = "marginal probability";
        rownames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[["caprese"]]$probabilities$probabilities.fit$estimated.conditional.probs) = "conditional probability";
                
    }
    else if(reconstruction$parameters$algorithm=="CAPRI") {
    	
    		cat("Executing now the estimation procedure, this may take a long time...\n")
    	
    		# go through the models
    	do.estimate.error.rates = FALSE;
    		if(is.na(error.rates[1])) {
    		do.estimate.error.rates = TRUE;
    	}
    		for (m in names(as.models(reconstruction))) {
    			
    			# if I also need to estimate the error rates
            if(do.estimate.error.rates) {
                # estimate the error rates
                error.rates = estimate.dag.error.rates(reconstruction$genotypes,as.marginal.probs(reconstruction,models=m)[[1]],as.joint.probs(reconstruction,models=m)[[1]],as.parents.pos(reconstruction,models=m)[[1]]);
            }
            
            # estimate the probabilities given the error rates
        		estimated.probabilities = estimate.dag.probs(reconstruction$genotypes,as.marginal.probs(reconstruction,models=m)[[1]],as.joint.probs(reconstruction,models=m)[[1]],as.parents.pos(reconstruction,models=m)[[1]],error.rates);
        	
        		# set the estimated error rates and probabilities
        		probabilities.fit = list(estimated.marginal.probs=estimated.probabilities$marginal.probs,estimated.joint.probs=estimated.probabilities$joint.probs,estimated.conditional.probs=estimated.probabilities$conditional.probs);
        
        		reconstruction$model[[m]]$error.rates = error.rates
        		reconstruction$model[[m]]$probabilities$probabilities.fit = probabilities.fit
        	
        		# set colnames and rownames
        		rownames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.marginal.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.marginal.probs) = "marginal probability";
            rownames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.joint.probs) = colnames(data$genotypes);
            rownames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.conditional.probs) = colnames(data$genotypes);
            colnames(reconstruction$model[[m]]$probabilities$probabilities.fit$estimated.conditional.probs) = "conditional probability";

    			
    		}
    }
    else {
    		stop("A valid algorithm has to be provided in order to estimate its confidence.",call.=FALSE);
    }
    
    reconstruction$parameters$do.estimation = TRUE;
    return(reconstruction);

}


#' Bootstrap a reconstructed progression model
#'
#' @examples
#' data(maf)
#' mutations = import.MAF(maf)
#' recon = tronco.capri(mutations)
#' boot = tronco.bootstrap(recon)
#' tronco.plot(boot)
#'
#' @title tronco bootstrap
#' @param reconstruction The output of tronco.capri or tronco.caprese
#' @param type TODO
#' @param nboot TODO
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.bootstrap
tronco.bootstrap <- function( reconstruction, 
                              type = "non-parametric", 
                              nboot = 100)
{
    # check for the inputs to be given
    if(is.null(reconstruction)) {
        stop("A valid reconstruction has to be provided in order to estimate its confidence.", call. = FALSE)
    }
    
    # check for the input to be compliant
    is.compliant(reconstruction)

    if(reconstruction$parameters$do.estimation == FALSE && type == "parametric") {
        stop("To perform parametric bootstrap, the estimation of the error rates and probabilities should be performed.", call. = FALSE)
    }

	if(type == "statistical" && !(reconstruction$parameters$algorithm == "CAPRI" && reconstruction$parameters$do.boot == TRUE)) {
        stop("To perform statistical bootstrap, the algorithm used for the reconstruction must by CAPRI with bootstrap.", call. = FALSE)
    }

    # set all the needed parameters to perform the bootstrap estimation
    if(type == "non-parametric" || type == "parametric" || type == "statistical") {
    	
        dataset = reconstruction$genotypes
        do.estimation = FALSE
        silent = TRUE
        
        if(!is.null(reconstruction$bootstrap)) {
        		bootstrap = reconstruction$bootstrap
        }
        else {
        		bootstrap = list()
        }
        
        if(reconstruction$parameters$algorithm == "CAPRESE") {
            lambda = reconstruction$parameters$lambda
        } else if(reconstruction$parameters$algorithm == "CAPRI") {
            
            if(!is.null(reconstruction$hypotheses)) {
                hypotheses = reconstruction$hypotheses
            } else {
            	hypotheses = NA
            }
            
            command.capri = reconstruction$parameters$command
            regularization = reconstruction$parameters$regularization
            do.boot = reconstruction$parameters$do.boot
            nboot.capri = reconstruction$parameters$nboot
            pvalue = reconstruction$parameters$pvalue
            min.boot = reconstruction$parameters$min.boot
            min.stat = reconstruction$parameters$min.stat
            boot.seed = reconstruction$parameters$boot.seed 
        }
    } else {
    	stop("The types of bootstrap that can be performed are: non-parametric, parametric or statistical.", call. = FALSE)
    }

    # perform the selected bootstrap procedure
    cat("Executing now the bootstrap procedure, this may take a long time...\n")
    cat("Expected completion in approx.",round((reconstruction$execution.time[3]*nboot)/60,digits=0),"minutes.\n")

    if(reconstruction$parameters$algorithm == "CAPRESE") {
        
        curr.boot = bootstrap.caprese(dataset,
<<<<<<< HEAD
        							  lambda,
        							  do.estimation,
                            		  silent,
                            		  reconstruction, 
                            		  type,
                            		  nboot)
=======
        							      lambda,
        							      do.estimation,
                            		      silent,
                            		      reconstruction, 
                            		      type,
                            		      nboot,
                            		      bootstrap)
>>>>>>> a212964d8bfb780b0a3abec02d48a98fcd2ca8e2
                                          
        reconstruction$bootstrap = curr.boot
        
        cat(paste("\nPerformed ", type, " bootstrap with ", nboot, " resampling and ", lambda, " as shrinkage parameter.\n\n", sep =""))
    
    } else if(reconstruction$parameters$algorithm == "CAPRI") {
    	
<<<<<<< HEAD
    	curr.boot = bootstrap.capri(dataset, 
                            		hypotheses, 
                            		command.capri, 
                            		regularization, 
                            		do.boot,
                            		nboot.capri, 
                            		pvalue,
                            		min.boot,
                            		min.stat,
                            		boot.seed,
                            		do.estimation,
                            		silent,
                            		reconstruction, 
                            		type,
                            		nboot)
=======
    		curr.boot = bootstrap.capri(dataset, 
                            		    hypotheses, 
                            		    command.capri, 
                            		    regularization, 
                            		    do.boot,
                            		    nboot.capri, 
                            		    pvalue,
                            		    min.boot,
                            		    min.stat,
                            		    boot.seed,
                            		    do.estimation,
                            		    silent,
                            		    reconstruction, 
                            		    type,
                            		    nboot,
                            			bootstrap)
>>>>>>> a212964d8bfb780b0a3abec02d48a98fcd2ca8e2

        reconstruction$bootstrap = curr.boot

        if(do.boot == TRUE) {
        	cat(paste("\nPerformed ", type, " bootstrap with ", nboot, " resampling and ", pvalue, " as pvalue for the statistical tests.\n\n", sep =""))
        } else {
        	cat(paste("\nPerformed ", type, " bootstrap with ", nboot, " resampling.\n\n", sep =""))
        }
    }

    return(reconstruction)
}

#### end of file -- tronco.bootstrap.R


#### tronco.plot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.




###########################
####### TRONCO PLOT #######
###########################



#' @export tronco.plot
#' @title plot a progression model
#'
#' @description
#' \code{tronco.plot} plots a progression model from a recostructed \code{curr.reconstruction}. 
tronco.plot = function(x,
                       regularization=names(x$model),
                       fontsize = NA, 
                       height=2,
                       width=3,
                       height.logic = 1,
                       pf = FALSE, 
                       disconnected=FALSE,
                       scale.nodes=NA,
                       name=deparse(substitute(capri)),
                       title = as.name(x),  
                       confidence = NA, 
                       p.min = x$parameters$pvalue,
                       legend = TRUE, 
                       legend.cex = 1.0, 
                       edge.cex = 1.0,
                       label.edge.size = NA, 
                       expand = TRUE,
                       genes = NULL,
                       relations.filter = NA,
                       edge.color = 'black',
                       pathways.color = 'Set1',
                       file = NA, # print to pdf,
                       legend.pos = 'bottom',
                       pathways = NULL,
                       lwd = 3,
                       ...
                       ) 
{
  hidden.and = F
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
    library(Rgraphviz)
  }

  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }
  
  # Checks if reconstruction exists
  if(missing(x)) {
    stop("reconstruction missing, usage: hypo.plot(reconstruction, ...", call.=FALSE);
  }
    
  logical_op = list("AND", "OR", "NOT", "XOR", "*", "UPAND", "UPOR", "UPXOR")
 

  if(length(regularization) > 2) {
    stop("Too many regularizators (max is 2)", call.=FALSE)
  }

  if(!regularization[1] %in% names(x$model)) {
    stop(paste(regularization[1], "not in model"), call.=FALSE);
  }

  #print(regularization)

  sec = ifelse(length(regularization) == 2, T, F)

  if(sec && !regularization[2] %in% names(x$model)) {
    stop(paste(regularization[2], "not in model"), call.=FALSE);
  }
  
  # Models objects
  primary = as.models(x, models = regularization[1])[[1]]    
  if(sec) secondary = as.models(x, models = regularization[2])[[1]]

  # USARE getters adj.matrix
  if (sec && !all( rownames(primary$adj.matrix$adj.matrix.fit) %in% rownames(secondary$adj.matrix$adj.matrix.fit))) {     
    stop("primary and secondary must have the same adj.matrix! See: the function tronco.bootstrap.", call.=FALSE)
  }
  
  # Get the adjacency matrix - this could have been donw with getters
  adj.matrix = primary$adj.matrix
  if(sec) adj.matrix = secondary$adj.matrix
  c_matrix = adj.matrix$adj.matrix.fit
    
  if(is.function(relations.filter))
  {
    cat('*** Filtering relations according to function "relations.filter", visualizing:\n')
    adj.matrix = as.adj.matrix(x, models = regularization)
    sel.relation = as.selective.advantage.relations(x, models = regularization)
    #sel.relation = lapply(sel.relation, keysToNames, x=x)
    
    # Select only relations which get TRUE by "relations.filter"
    sel.relation = lapply(sel.relation, 
                          function(z){
                            # apply can not be used - implicit cohersion to char is crap
                            # z[ apply(z, 1, relations.filter), ]
                            ####
                            mask = rep(T, nrow(z))
                            for(i in 1:nrow(z))
                              mask[i] = relations.filter(z[i, ]) 
                            return(z[mask, , drop = F])                              
                            })
    
    print(sel.relation)
    
    sel.relation = get(regularization[2], sel.relation)
    
    
    c_matrix.names = rownames(c_matrix)
    c_matrix = matrix(0, nrow = nrow(c_matrix), ncol = ncol(c_matrix))
    rownames(c_matrix) = c_matrix.names
    colnames(c_matrix) = c_matrix.names
    
    cat(paste0('Selected relations: ', nrow(sel.relation), ' [out of ', nrow(as.selective.advantage.relations(x, models = regularization)[[2]]), ']\n'))
    
    if(nrow(sel.relation) > 0)
    for(i in 1:nrow(sel.relation))
      c_matrix[ nameToKey(x, sel.relation[i, 'SELECTS']), nameToKey(x, sel.relation[i, 'SELECTED'])] = 1
      
    #print(lapply(adj.matrix, keysToNames, x=x))
    
  }
  
  # get the probabilities
  probabilities = primary$probabilities
  if(sec) {
    probabilities = secondary$probabilities
  }
  marginal_p = probabilities$probabilities.observed$marginal.probs
  
  # if prima facie change the adj matrix
  if (pf) {
    c_matrix = adj.matrix$adj.matrix.pf
  }
  
  if (all(c_matrix == F) || (sec && all(primary$adj.matrix$adj.matrix.fit == F))) {
    warning('No edge in adjacency matrix! Nothing to show here.')
    return(NULL)
  }
  
  #bootstrap = x$bootstrap
  #if(sec) {
  #  bootstrap = secondary$bootstrap
  #}

  #if (confidence) {
  #  conf_matrix = if (pf) bootstrap$edge.confidence$edge.confidence.pf else bootstrap$edge.confidence$edge.confidence.bic
  #}
  #print(c_matrix)
  # conf_matrix = NULL
  
  # get algorithm parameters
  parameters = x$parameters
  # print(parameters)
  
  # get hypotheses
  hypotheses = x$hypotheses
  hstruct = NULL
  if (!is.null(hypotheses) && !is.na(hypotheses) ) {
    hstruct = hypotheses$hstructure
  }

  # get event from genes list
  events = NULL
  if (is.vector(genes)) {
    events = unlist(lapply(genes, function(x){names(which(as.events(x)[,'event'] == x))}))
  }
  
  cat('*** Expanding hypotheses syntax as graph nodes:')

  # expand hypotheses
  #if (is.na(confidence)) {
    expansion = hypotheses.expansion(c_matrix, 
                                     hstruct, 
                                     hidden.and, 
                                     expand, 
                                     events)
    hypo_mat = expansion[[1]]
    hypos_new_name = expansion[[2]]
  #} else {
  #  expansion = hypotheses.expansion(c_matrix, 
  #                                   hstruct, 
  #                                   hidden.and, 
  #                                   expand, 
  #                                   events, 
  #                                   conf_matrix)
  #  hypo_mat = expansion[[1]]
  #  hypos_new_name = expansion[[2]]
  #  conf_matrix = expansion[[3]]
  #}
  
  # print(hypo_mat)
  
  cat('\n*** Rendering graphics\n')
  
  # remove disconnected nodes
  if(!disconnected) { 
    cat('Nodes with no incoming/outgoing edges will not be displayed.\n')
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
  }
  
  
  attrs = list(node = list())
      
  #print(hypo_mat)
  
  hypo_graph = graph.adjacency(hypo_mat)
  #cat('\n')
  #print(V(hypo_graph)$name[26:30])
  v_names = gsub("_.*$", "", V(hypo_graph)$name)
  if (!expand) {
    v_names = gsub("^[*]_(.+)", "*", V(hypo_graph)$name)
  }
  #print(v_names[26:30])
  new_name = list()
  
  for(v in v_names) {
    if(v %in% rownames(x$annotations)) {
      n = x$annotations[v,"event"]
      new_name = append(new_name, n)
    } else {
      new_name = append(new_name, v)
    }
  }
    
  #print(V(hypo_graph)$name)
  #print(v_names)
  #print(new_name)
  #print(V(hypo_graph)$label)
  #print(hypo_graph)
  
  
  V(hypo_graph)$label = new_name
  graph <- igraph.to.graphNEL(hypo_graph)
  
  node_names = graph::nodes(graph)
  nAttrs = list()
  
  nAttrs$label = V(hypo_graph)$label
  names(nAttrs$label) = node_names
  
  # set a default color
  nAttrs$fillcolor =  rep('White', length(node_names))
  names(nAttrs$fillcolor) = node_names
  
  # set fontsize
  
  if(is.na(fontsize)) {
    fontsize = 24 - 4*log(nrow(hypo_mat))
    cat(paste0('Set automatic fontsize scaling for node labels: ', fontsize, '\n'))
  }
  nAttrs$fontsize = rep(fontsize, length(node_names))
  names(nAttrs$fontsize) = node_names

  # set node shape
  nAttrs$shape = rep('ellipse', length(node_names))
  names(nAttrs$shape) = node_names

  # set node height
  nAttrs$height = rep(height, length(node_names))
  names(nAttrs$height) = node_names
  
  # set node width
  nAttrs$width = rep(width, length(node_names))
  names(nAttrs$width) = node_names


  
  short.label = nAttrs$label
  names(short.label) = names(nAttrs$label)
  if (!is.na(scale.nodes)) {

    
    # foreach node
    min_p = min(marginal_p)
    max_p = max(marginal_p)
    #print(scale.nodes)

    for (node in node_names) {
      prefix = gsub("_.*$", "", node)
      if ( !(prefix %in% logical_op)) {
        # Scaling ANDRE
         increase_coeff = scale.nodes + (marginal_p[node,] - min_p) / (max_p - min_p)
         nAttrs$width[node] = nAttrs$width[node] * increase_coeff
         nAttrs$height[node] = nAttrs$height[node] * increase_coeff
         nAttrs$label[node] = paste0(nAttrs$label[node], '\\\n', round(marginal_p[node, ]*100, 0), '%', ' (', sum(as.genotypes(x)[, node]) ,')')

        # increase_coeff = (1- (max_p - marginal_p[node,])) * scale.nodes
        
# #         print('***')
        # print(marginal_p[node,])
        # print(max_p)
        # print(increase_coeff)
        
        # nAttrs$width[node] = nAttrs$width[node] + increase_coeff
        # nAttrs$height[node] = nAttrs$height[node] + increase_coeff
        
        # print(nAttrs$width[node])        
      }
    }
  }
  
  # use colors defined in tronco$types
    w = unlist(lapply(names(nAttrs$fillcolor), function(w){
    if (w %in% rownames(x$annotations)) {
      x$types[x$annotations[w,'type'], 'color']
     }
    else
      'White'
    }))
  nAttrs$fillcolor[] = w

  legend_logic = NULL
  
  # set color, size form and shape each logic nodes (if hypos expansion actived)
  node.type = 'box'
  if (expand) {
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
    if (any(w)) {
      legend_logic['Exclusivity (soft)'] = 'orange'
    }
    nAttrs$fillcolor[which(w)] = 'orange'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
    if (any(w)) {
      legend_logic['Co-occurence'] = 'darkgreen'
    }
    nAttrs$fillcolor[which(w)] = 'darkgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
    if (any(w)) {
      legend_logic['Exclusivity (hard)'] = 'red'
    }
    nAttrs$fillcolor[which(w)] = 'red'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic

    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPOR'
    if (any(w)) {
      legend_logic['Exclusivity (soft)'] = 'orange'
    }
    nAttrs$fillcolor[which(w)] = 'orange'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPAND'
    if (any(w)) {
      legend_logic['Co-occurence'] = 'lightgreen'
    }
    nAttrs$fillcolor[which(w)] = 'lightgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPXOR'
    if (any(w)) {
      legend_logic['Exclusivity (hard)'] = 'red'
    }
    nAttrs$fillcolor[which(w)] = 'red'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic

  }
  #print(legend_logic)
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == '*'
  if (any(w)) {
      legend_logic['Co-occurence'] = 'darkgreen'
    }
  nAttrs$fillcolor[which(w)] = 'darkgreen'
  nAttrs$label[which(w)] = ''
  nAttrs$shape[which(w)] = node.type
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = height.logic

  
  # node border thickness
  #nAttrs$lwd = 1

  # node border to black
  nAttrs$color = rep("black", length(node_names))
  names(nAttrs$color) = node_names

  nAttrs$fontcolor = rep("black", length(node_names))
  names(nAttrs$fontcolor) = node_names

  # set node border based on pathways information
  #cat('\n')
  legend_pathways = NULL
  if(!is.null(pathways)) {
    cat('Annotating nodes with pathway information. \n')
    
    if(length(pathways.color) == 1 && pathways.color %in% rownames(brewer.pal.info)) 
    {
        cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
        cols = brewer.pal(n=length(names(pathways)), name=pathways.color)
    }
    else
    {
        if(length(pathways.color) != length(names(pathways))) 
            stop('You did not provide enough colors to annotate', length(names(pathways)), 'pathways. 
                    Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')
        cols = pathways.color
    }
    
    # print(cols)   
    # print(names(cols))
    # print(names(pathways))
        
    
    names(cols) = names(pathways)
    #print(cols)

    #nAttrs$lwd = lwd
    #nAttrs$col = rep("white", length(node_names))
    names(nAttrs$col) = node_names

    for(path in names(pathways)) {
      #cat('\npath: ', pathways[[path]])
      n = short.label[which(short.label %in% pathways[[path]])]
      #cat('\nfound: ', unlist(n), unlist(names(n)))
      #cat('color: ', cols[[path]])
      #cat('\n')
      nAttrs$color[unlist(names(n))] = cols[[path]]
      nAttrs$fontcolor[unlist(names(n))] = cols[[path]]
      
      
      if(length(n) > 0) {
        legend_pathways[path] = cols[[path]]
      }
    }
  }

  
  # edges properties
  
  edge_names = edgeNames(graph)
  eAttrs = list()
  
  # set temporary edge shape
  eAttrs$lty = rep("solid", length(edge_names))
  names(eAttrs$lty) = edge_names
  
  #set edge thikness based on prob
  eAttrs$lwd = rep(1, length(edge_names))
  names(eAttrs$lwd) = edge_names
  
  #set edge name based on prob
  eAttrs$label = rep('', length(edge_names))
  names(eAttrs$label) = edge_names
  
  #set fontsize to label.edge.size (default)
  if(is.na(label.edge.size)) {
    label.edge.size = fontsize/2      
    cat(paste0('Set automatic fontsize for edge labels: ', label.edge.size, '\n'))    
  }
  eAttrs$fontsize = rep(label.edge.size, length(edge_names))
  names(eAttrs$fontsize) = edge_names
  
  #set edge color to black (default)
  eAttrs$color = rep(ifelse(sec, 'darkgrey', edge.color), length(edge_names))
  names(eAttrs$color) = edge_names

  #set edge arrowsize to 1 (default)
  eAttrs$arrowsize = rep(1 * edge.cex, length(edge_names))
  names(eAttrs$arrowsize) = edge_names
  
  #record logic edge
  eAttrs$logic = rep(F, length(edge_names))
  names(eAttrs$logic) = edge_names
  

  #print('confidence')
  #print(confidence)
  
  if(any(!is.na(confidence))) {
    cat('Adding confidence information: ')
    conf = as.confidence(x, confidence)
    
    cat(paste(paste(confidence, collapse = ', '), '\n'))
    # names.conf = names(conf)
    
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
    
      from = edge[1]
      to = edge[2]

      #cat("***\n", from, '->', to, '\n')
      
      pval.names = c('hg', 'pf', 'tp')
      red.lable = FALSE

      if(is.logic.node.up(from) || is.logic.node.down(to)) {
        next
      }

      if(from %in% names(hypos_new_name)){ conf_from = hypos_new_name[[from]] } else { conf_from = from }
      if(to %in% names(hypos_new_name)){ conf_to = hypos_new_name[[to]] } else { conf_to = to }



        #cat(conf_from, '->', conf_to, '\n')
      
      for(i in confidence) {
          
        conf_p = get(i, as.confidence(x, i))
        if(! (conf_from %in% rownames(conf_p) && conf_to %in% colnames(conf_p))) {
          next
        }
          
        eAttrs$label[e] = paste0(
        eAttrs$label[e],
        ifelse(conf_p[conf_from, conf_to] < 0.01, "< 0.01", round(conf_p[conf_from, conf_to], 2)))
        #conf_p[conf_from, conf_to], '\\\n')
        #cat(conf_p[conf_from, conf_to], '\n')
        
        if(i %in% pval.names && conf_p[conf_from, conf_to] > p.min) {
          eAttrs$fontcolor[e] = 'red'
          eAttrs$label[e] = paste0(eAttrs$label[e], ' *')
        }
        eAttrs$label[e] = paste0(eAttrs$label[e], '\\\n')
                  
      }
      
      
# #    
      # # ..checks if confidence is available
      # if (from %in% rownames(conf_matrix) && to %in% colnames(conf_matrix)) {
        # if(from %in% names(hypos_new_name)){ conf_from = hypos_new_name[[from]] } else { conf_from = from }
        # if(to %in% names(hypos_new_name)){ conf_to = hypos_new_name[[to]] } else { conf_to = to }
        # # if confidence > 0..
        # # print(paste('from', from, 'to', to, ':', conf_matrix[from, to]))
        
        
        # if (conf_matrix[from, to] == 1) {
          # # ..set edge thickness and label..
          # eAttrs$label[e] = '100%'
          # eAttrs$lwd[e] = log(150)
        # } else if (conf_matrix[from, to] >= 0.01) {
          # # ..draw it on the graph..

          # eAttrs$label[e] = paste0('', round(conf_matrix[from, to] * 100, 0), '%')
          # eAttrs$lwd[e] = log(conf_matrix[from, to] * 150)
        # } else {
          # # ..else set the style of the edge to dashed
          # eAttrs$label[e] = "< 1%"
          # eAttrs$lwd[e] = log(1.5)
        # }


    #print(conf_from)
    #print(conf_to)
    
    
        # hyper_geom = x$confidence[[3]][conf_from, conf_to]
        # if (hyper_geom < 0.01) { hyper_geom = '< .01'} else { hyper_geom = round(hyper_geom, 2)}
        # eAttrs$label[e] = paste(eAttrs$label[e], '  ', hyper_geom)


      # } else {
        # # ..else this edge is located inside to an hypothesis, so no arrow to show
        # eAttrs$logic[e] = T
        # eAttrs$color[e] = 'black'
      # }
    }
    cat('RGraphviz object prepared.\n')
  }

  # remove arrows from logic node (hidden and)
  for(e in edge_names) {
    edge = unlist(strsplit(e, '~'))
    from = edge[1]
    to = edge[2]
    
    # print(e)
    # print(is.logic.node.down(to))
    # print(is.logic.node.up(from))
        
    if (is.logic.node.down(to)) {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
     
    # print(from)
    # print(substr(to, start=1, stop=2))
     # eAttrs$color[e] = 'red'
     
     if(substr(to, start=1, stop=2) == 'OR')
         eAttrs$color[e] = 'orange'
     if(substr(to, start=1, stop=3) == 'XOR')
         eAttrs$color[e] = 'red'
     if(substr(to, start=1, stop=3) == 'AND')
         eAttrs$color[e] = 'darkgreen'
           
     eAttrs$lty[e] = 'dashed'

    nAttrs$shape[to] = 'circle' 
    # nAttrs$shape[to] = function(w, ...) {
        # for(i in 1:3) points(w, cex=i, ...)
    # }
    
    # print(nAttrs)

    }

    if (is.logic.node.up(from)) {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
      # eAttrs$color[e] = 'black'
      
      eAttrs$lty[e] = 'dashed'
      
        # print(substr(to, start=1, stop=2))
        # print(substr(from, start=1, stop=2))
    # print(from)
    # print(to)
 
 
      if(substr(from, start=1, stop=4) == 'UPOR')
         eAttrs$color[e] = 'orange'
     if(substr(from, start=1, stop=5) == 'UPXOR')
         eAttrs$color[e] = 'red'
     if(substr(from, start=1, stop=5) == 'UPAND')
         eAttrs$color[e] = 'darkgreen'
     
      
    } else if(substr(from, start=1, stop=1) == '*') {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
      eAttrs$color[e] = 'black'
    }
  }
  
  #print(eAttrs$lty)
  
  if(pf) {
    cat('*** Add prima facie edges: ')
    # for each edge..
    bic = adj.matrix$adj.matrix.bic
    #print(bic)
    #print('logic edge')
    #print(eAttrs$logic)
    
    #print(rownames(bic))
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      old_name = hypos_new_name[[from]]
      #cat('\n\nnodo from:', hypos_new_name[[from]])
      if (!is.null(old_name)) {
        from = old_name
      }
      to = edge[2]
      if (substr(to, start=1, stop=1) == '*') {
        to = substr(to, start=3, stop=nchar(to))
      }
      
      # ..checks if edge is present in BIC

      #cat('\nfrom:', from, 'to:', to)
      #cat('\nfrom in bic? ', (from %in% rownames(bic)))
      #cat('\nto in bic? ', (to %in% rownames(bic)))
      #cat('\n!eAttrs$logic[e]', !eAttrs$logic[e])
      # check if edge in BIC (valid only if not logic edge) and 'to' is not a fake and
      if ( (from %in% rownames(bic)) &&
           (to %in% colnames(bic)) &&
           !eAttrs$logic[e] &&
           bic[from, to] == 0
           ) {
        #cat("\nprima facie!!!")
        eAttrs$color[e] = 'red'
      } else {
        #cat('\nno PF!')
      }
    }
    cat('done')
  }

  if (sec) {
    #print(edge_names)
    
    pri.adj = primary$adj.matrix$adj.matrix.fit
    #print(rownames(pri.adj))
    #print(hypos_new_name)
    for(from in rownames(pri.adj)) {
      for(to in colnames(pri.adj)) {
        from.alt.name = from
        if (from %in% hypos_new_name) {
          from.alt.name = names(which(hypos_new_name == from))
        }

        if(pri.adj[from, to] == 1) {
          eAttrs$color[paste(from.alt.name, to, sep='~')] = edge.color
        }


        #cat('\nfrom: ', from.alt.name, ' to: ', to)
      }
    }
  }
  

  
  # tp = x$confidence['temporal priority', ][[1]]
  # pr = x$confidence['probability raising', ][[1]] 
  # hg = x$confidence['hypergeometric test', ][[1]] 
  # # conf = apply(tp, pr, hg)
  
  
  
  # # print(eAttrs$label)
  
  # p.min = 0.01
  
  # for(i in names(eAttrs$label))
  # {
    # names = strsplit(i, '~')[[1]]
    # print(i)
    # print(names)
    
    # if(all(names %in% colnames(tp)))      
    # {
      # val = 1 - round(max(
        # # tp[names[1], names[2]]
        # # ,
        # # pr[names[1], names[2]],
        # hg[names[1], names[2]]
        # )
        # , 3)
        
        # print(val)
      
      # val = ifelse(val >= (1 - p.min), 1, val)  
      
      # if(val < 1)   
      # {
        # eAttrs$label[i] = paste(val, '\\\n', 'sss')
        # eAttrs$fontcolor[i] = 'red'
        
      # # eAttrs$lwd[i] = eAttrs$lwd[i] * exp(val) * 1.5
      # # eAttrs$color[i] = 'red'
      # }
    # }
    
    
    
  # }
      #print(eAttrs)

  
  
  # print(eAttrs)
  
  cat('Plotting graph and adding legends.\n')
  plot(graph, nodeAttrs=nAttrs, edgeAttrs=eAttrs, main=title, ... )
  
  # Adds the legend to the plot
  if (legend) {
    valid_events = colnames(hypo_mat)[which(colnames(hypo_mat) %in% colnames(c_matrix))]
    legend_names = unique(x$annotations[which(rownames(x$annotations) %in% valid_events), 'type'])
    pt_bg = x$types[legend_names, 'color']
    legend_colors = rep('black', length(legend_names))
    pch = rep(21, length(legend_names))
    
    if (length(legend_logic) > 0) {
      pch = c(pch, 0, 0, rep(22, length(legend_logic)))
      legend_names = c(legend_names, ' ', expression(bold('Patterns')), names(legend_logic))
      legend_colors = c(legend_colors, 'white', 'white', rep('black', length(legend_logic)))
      pt_bg = c(pt_bg, 'white', 'white', legend_logic)  
    }

    if (length(legend_pathways) > 0) {
      pch = c(pch, 0, 0, rep(21, length(legend_pathways)))
      legend_names = c(legend_names, ' ', expression(bold('Pathways')), names(legend_pathways))
      pt_bg = c(pt_bg, 'white', 'white', rep('white', length(legend_pathways)))
      legend_colors = c(legend_colors, 'white', 'white', legend_pathways)  
    }
    
    if (legend.pos == 'bottom') {
      legend.pos.l = 'bottomleft'
      legend.pos.r = 'bottomright'
    } else if (legend.pos == 'top') {
      legend.pos.l = 'topleft'
      legend.pos.r = 'topright'
    } else {
      legend.pos.l = locator(1)
      legend.pos.r = locator(1)
    }

    legend(legend.pos.r,
           legend = legend_names,
           title = expression(bold('Events type')),
           bty = 'n',
           cex = legend.cex,
           pt.cex = 1.5 * legend.cex,
           pch = pch,
           col = legend_colors,
           pt.bg = pt_bg)

    #add thickness legend
    valid_names = node_names

    if(!disconnected) { 
      del = which(rowSums(hypo_mat) + colSums(hypo_mat) == 0 )
      w = !(rownames(hypo_mat) %in% names(del))
      valid_names = rownames(hypo_mat[w,])
    }

    if(expand) {
      valid_names = valid_names[unlist(lapply(valid_names, function(x){!is.logic.node(x)}))]
    }

    valid_names = grep('^[*]_(.+)$', valid_names, value = T, invert=T)

    #dim = nAttrs$height[valid_names]
    #prob = marginal_p[valid_names, ]
    
    #min = min(dim)
    #p_min = round(min(prob) * 100, 0)
    #max = max(dim)
    #p_max = round(max(prob) * 100, 0)
    
    
    #marginal_p = marginal_p[valid_names, , drop = FALSE]

    # Get label of the (first) event with minimum marginale 
    #min.p =   rownames(marginal_p)[which(min(marginal_p) == marginal_p) ]
    #label.min = as.events(x)[ min.p[1] , , drop = FALSE]

    # Get label of the (first) event with max marginale 
    #max.p = rownames(marginal_p)[which(max(marginal_p) == marginal_p) ]
    #label.max = as.events(x)[ max.p[1] , ,  drop = FALSE]

    # Frequency labels
    #min.freq = round(min(marginal_p) * 100, 0)
    #max.freq = round(max(marginal_p) * 100, 0)
      
    #freq.labels = c( 
      #paste0(min.freq, ifelse((min.freq < 10 && max.freq > 9), '%  ', '%'), ' ', label.min[, 'event'], ' (min)'),
      #paste0(max.freq, '% ', label.max[, 'event'], ' (max)')
    #)

    freq.labels = ""
    stat.pch = 0
    pt.bg = "white"
    col = "white"
    if (any(!is.na(confidence))) {
      freq.labels = c(expression(bold('Edge confidence')), lapply(confidence, function(x){
        if(x == "hg")
          return("Hypergeometric test")
        if(x == "tp")
          return("Temporal Priority")
        if(x == "pr")
          return("Probability Raising")
        if(x == "pb")
          return("Parametric Bootstrap")
        if(x == "npb")
          return("Non Parametric Bootstrap")

      }), paste("p <", p.min))
      #stat.pch = c(21, 21)
      stat.pch = c(0, rep(18, length(confidence)), 0)
      pt.bg = c('white', rep('white', length(confidence)), 'white')
      col = c('white', rep('black', length(confidence)), 'white')
    }
  



        
  if('Pattern' %in% as.types(x)) 
            y = delete.type(x, 'Pattern')
    else y = x

            
    freq.labels = c(freq.labels, 
      ' ',
      expression(bold('Sample size')),
      paste0('n = ', nsamples(x), ', m = ', nevents(x)),
      paste0('|G| = ', ngenes(y), ', |P| = ', npatterns(x))     
    ) 
    
    reg.labels = c( '\n',
      expression(bold('Regularization')),
      paste0(names(x$model))
    )


     
    stat.pch = c(stat.pch, rep(0, 2), rep(20, 2), rep(0, 2), rep(20, 2))
    pt.bg = c(pt.bg, rep('white', 2), rep('black', 2), rep('white', 2), 'black', 'darkgrey')
    col = c(col, rep('white', 2), rep('white', 2), rep('white', 2),'black', 'darkgrey') 
    
    # print(data.frame(stat.pch, pt.bg, col))
    
    legend(legend.pos.l,
           legend = c(freq.labels, reg.labels),
           title = "",
           bty = 'n',
           box.lty = 3,
           box.lwd = .3,
           pch = stat.pch,
           pt.cex = 1.5  * legend.cex,
           ncol = 1,
           pt.bg = pt.bg,
           cex = legend.cex,
           col = col)
       
  }

  
  if(!is.na(file))
  {
    dev.copy2pdf(file = file)
  }
  cat('\n')
}

#' @import Rgraphviz
#' @import igraph
#' @import RColorBrewer
tronco.consensus.plot = function(models,
                       secondary=NULL, 
                       fontsize=18, 
                       MIN.HITS = 1,
                       height=2,
                       width=3,
                       height.logic = 1,
                       pf = FALSE, 
                       disconnected=FALSE,
                       scale.nodes=NA,
                       title = paste("Consensus Progression model"),  
                       confidence = FALSE, 
                       legend = TRUE, 
                       legend.cex = 1.0, 
                       edge.cex = 1.0,
                       label.edge.size = 12, 
                       hidden.and = T,
                       expand = TRUE,
                       edge.color = 'black',
                       pathways.color = 'Set1',
                       file = NA, # print to pdf,
                       legend.pos = 'bottom',
                       pathways = NULL,
                       ...
                       ) 
{
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
    library(Rgraphviz)
  }

  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
  }
  
  # Checks if reconstruction exists
  if(missing(models) || !is.list(models)) {
    stop("Models missing, usage: ... ...", call.=FALSE);
  }
    
  logical_op = list("AND", "OR", "NOT", "XOR", "*", "UPAND", "UPOR", "UPXOR")
  
  
smaller.to.bigger = function(m,cn)
{
    x = matrix(0, nrow = length(cn), ncol = length(cn))
    rownames(x) = cn
    colnames(x) = cn
 
    
    for(i in 1:nrow(m))
        for(j in 1:nrow(m))
            x[rownames(m)[i], rownames(m)[j]] = ifelse(m[i,j] == 1, 1, 0) 
    return(x)
}


    # All the adjacency matrices
    matrices = list()
    for(i in 1:length(models))
        matrices = append(matrices, list(models[[i]]$adj.matrix$adj.matrix.bic))

    # All their colnames - all possible eventsand types
    cn = unique(Reduce(union, lapply(matrices, colnames)))
    
    all.events = NULL
    for(i in 1:length(models)) all.events = rbind(all.events, as.events(models[[i]]$data))
    all.events = unique(all.events)

    all.types = NULL
    for(i in 1:length(models)) all.types = rbind(all.types, models[[i]]$data$types)
    all.types = unique(all.types)
    
    # Consensus + overall adjacency matrix
    consensus = Reduce('+', lapply(matrices, smaller.to.bigger, cn=cn))
    adjacency = consensus
    adjacency[adjacency < MIN.HITS] = 0 
    adjacency[adjacency > 1] = 1 
        
    cat('Consensus adjacency matrix:', nrow(adjacency), 'x', ncol(adjacency), ', minimum consensus', MIN.HITS, '\n')

    # All the marginal probabilities - entries from each input model
    marginal.probabilities = matrix(nrow = nrow(adjacency), ncol = 1)
    rownames(marginal.probabilities) = rownames(adjacency)
  
    for(i in 1:length(models))
    {
        model.marginals = models[[i]]$probabilities$probabilities.bic$marginal.probs
        marginal.probabilities[rownames(model.marginals),] = model.marginals
    }

    # Algorithm parameters should be the same for each model
    parameters = models[[1]]$parameters
  
    # Retrieve all hypotheses - TODO constraint this to a model with hypotheses
    all.hypos = list()
    for(i in 1:length(models))
        all.hypos = append(all.hypos, models[[i]]$data$hypotheses$hstructure)
    
    # We unpack the env and collapse it back
    all.hypos = sapply(all.hypos, as.list) # make them as list
    all.hypos.names = NULL
    hstruct = list()
    for(i in 1:length(all.hypos))
    {
        # Get hypotheses missing in hstruct
        local.env = all.hypos[[i]]
        local.env.names = names(local.env)[which( !(names(local.env) %in% names(hstruct)) )]
        
        all.hypos.names = c(all.hypos.names, local.env.names)
        hstruct = append(hstruct, local.env[local.env.names])
    }
    
    cat('Found', length(hstruct), 'hypotheses\n')
    hstruct = as.environment(hstruct)

    # print('originale ordine') 
         # print(colnames(adjacency))
     # print('**')
    # print(colnames(consensus))
     # print('**')


    # hypotheses.expansion requires hypotheses to be stored as rightmost columns
    adjacency = adjacency[, c(setdiff(colnames(adjacency), all.hypos.names), all.hypos.names)]
    adjacency = adjacency[c(setdiff(colnames(adjacency), all.hypos.names), all.hypos.names), ]
    
    consensus = consensus[, c(setdiff(colnames(consensus), all.hypos.names), all.hypos.names)]
    consensus = consensus[c(setdiff(colnames(consensus), all.hypos.names), all.hypos.names), ]
    
    all.events = all.events[colnames(adjacency), ]

    # print('nuovo ordine') 

         # print(colnames(adjacency))
     # print('**')
    # print(colnames(consensus))
     # print('**')
    
    # Expand hypotheses
    expansion = hypotheses.expansion(adjacency, 
                                     map = hstruct, 
                                     hidden_and = hidden.and, 
                                     expand = expand)
    hypo_mat = expansion[[1]]
    hypos_new_name = expansion[[2]]
    # print(hypos_new_name)
    
    # print('** espansa')
    # print(colnames(hypo_mat))

    # print('\n\n\n\n')

    # df = data.frame(rep('NA', ncol(adjacency)), row.names=paste(1:ncol(adjacency)))
    # df$adj = colnames(adjacency)
    # df$con = colnames(consensus)
    # df$hypo = colnames(hypo_mat)
    
    # print(df)

  # Remove disconnected nodes
  if(!disconnected) { 
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
  }
  
  cat('\n*** Render graphics: ')
  
  attrs = list(node = list())
      
  hypo_graph = graph.adjacency(hypo_mat)

  v_names = gsub("_.*$", "", V(hypo_graph)$name)
  if (!expand) {
    v_names = gsub("^[*]_(.+)", "*", V(hypo_graph)$name)
  }

  new_name = list()
  for(v in v_names) {
    if(v %in% rownames(all.events)) {
      n = all.events[v,"event"]
      new_name = append(new_name, n)
    } else {
      new_name = append(new_name, v)
    }
  }

  V(hypo_graph)$label = new_name
  graph <- igraph.to.graphNEL(hypo_graph)
  
  node_names = graph::nodes(graph)
  nAttrs = list()
  
  nAttrs$label = V(hypo_graph)$label
  names(nAttrs$label) = node_names
  
  # set a default color
  nAttrs$fillcolor =  rep('White', length(node_names))
  names(nAttrs$fillcolor) = node_names
  
  # set fontsize
  nAttrs$fontsize = rep(fontsize, length(node_names))
  names(nAttrs$fontsize) = node_names

  # set node shape
  nAttrs$shape = rep('ellipse', length(node_names))
  names(nAttrs$shape) = node_names

  # set node height
  nAttrs$height = rep(height, length(node_names))
  names(nAttrs$height) = node_names
  
  # set node width
  nAttrs$width = rep(width, length(node_names))
  names(nAttrs$width) = node_names
  
  if (!is.na(scale.nodes)) {

    min_p = min(marginal.probabilities)
    max_p = max(marginal.probabilities)

    for (node in node_names) {
      prefix = gsub("_.*$", "", node)
      if ( !(prefix %in% logical_op)) {
        # Scaling ANDRE
         increase_coeff = scale.nodes + (marginal.probabilities[node,] - min_p) / (max_p - min_p)
         nAttrs$width[node] = nAttrs$width[node] * increase_coeff
         nAttrs$height[node] = nAttrs$height[node] * increase_coeff    
      }
    }
  }
    
  # use colors defined in tronco$types
  w = unlist(lapply(names(nAttrs$fillcolor), function(x){
    if (x %in% rownames(all.events))
      all.types[all.events[x,'type'], 'color']
    else
      'White'
    }))
  nAttrs$fillcolor[] = w
  

  legend_logic = NULL
  
  # set color, size form and shape each logic nodes (if hypos expansion actived)
  node.type = 'box'
  if (expand) {
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
    if (any(w)) {
      legend_logic['Exclusivity (soft)'] = 'orange'
    }
    nAttrs$fillcolor[which(w)] = 'orange'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
    if (any(w)) {
      legend_logic['Co-occurence'] = 'darkgreen'
    }
    nAttrs$fillcolor[which(w)] = 'darkgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
    if (any(w)) {
      legend_logic['Exclusivity (hard)'] = 'red'
    }
    nAttrs$fillcolor[which(w)] = 'red'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic

    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPOR'
    if (any(w)) {
      legend_logic['Exclusivity (soft)'] = 'orange'
    }
    nAttrs$fillcolor[which(w)] = 'orange'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPAND'
    if (any(w)) {
      legend_logic['Co-occurence'] = 'lightgreen'
    }
    nAttrs$fillcolor[which(w)] = 'lightgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'UPXOR'
    if (any(w)) {
      legend_logic['Exclusivity (hard)'] = 'red'
    }
    nAttrs$fillcolor[which(w)] = 'red'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic

  }
  #print(legend_logic)
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == '*'
  if (any(w)) {
      legend_logic['Co-occurence'] = 'lightgreen'
    }
  nAttrs$fillcolor[which(w)] = 'lightgreen'
  nAttrs$label[which(w)] = ''
  nAttrs$shape[which(w)] = node.type
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = height.logic
  
  # node border to black
  nAttrs$color = rep("black", length(node_names))
  names(nAttrs$color) = node_names

  nAttrs$fontcolor = rep("black", length(node_names))
  names(nAttrs$fontcolor) = node_names

  # set node border based on pathways information
  legend_pathways = NULL
  if(!is.null(pathways)) {
    
    
    if(length(pathways.color) == 1 && pathways.color %in% rownames(brewer.pal.info)) 
    {
        cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
        cols = brewer.pal(n=length(names(pathways)), name=pathways.color)
    }
    else
    {
        if(length(pathways.color) != length(names(pathways))) 
            stop('You did not provide enough colors to annotate', length(names(pathways)), 'pathways. 
                    Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')
        cols = pathways.color
    }
    
    names(cols) = names(pathways)
    names(nAttrs$col) = node_names


    for(path in names(pathways)) {
      n = nAttrs$label[which(nAttrs$label %in% pathways[[path]])]
      nAttrs$color[unlist(names(n))] = cols[[path]]
      nAttrs$fontcolor[unlist(names(n))] = cols[[path]]
      if(length(n) > 0) {
        legend_pathways[path] = cols[[path]]
      }
    }

  }
  
  # edges properties
  
  edge_names = edgeNames(graph)
  eAttrs = list()
  
  # print(edge_names)
  
  # set temporary edge shape
  eAttrs$lty = rep("solid", length(edge_names))
  names(eAttrs$lty) = edge_names
  
  #set edge thikness based on prob
  eAttrs$lwd = rep(1, length(edge_names))
  names(eAttrs$lwd) = edge_names
  
  #set edge name based on prob
  eAttrs$label = rep('', length(edge_names))
  names(eAttrs$label) = edge_names
  
  #set fontsize to label.edge.size (default)
  eAttrs$fontsize = rep(label.edge.size, length(edge_names))
  names(eAttrs$fontsize) = edge_names
  
  #set edge color to black (default)
  eAttrs$color = rep(edge.color, length(edge_names))
  names(eAttrs$color) = edge_names

  #set edge arrowsize to 1 (default)
  eAttrs$arrowsize = rep(1 * edge.cex, length(edge_names))
  names(eAttrs$arrowsize) = edge_names
  
  #record logic edge
  eAttrs$logic = rep(F, length(edge_names))
  names(eAttrs$logic) = edge_names
  
  # cat('done')
   
  
  # print('RINOCERONTE') 
   
    # for each edge..
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      to = edge[2]
      # print('*(***)')
      # cat('=====FROM', from, '\n')
      # cat('to', to, '\n')
      # print(consensus[from, to])
      # print(sort(colnames(consensus)))
      # print(rownames(consensus))
      # idx.from = which(rownames())
      
      
      if(from %in% names(hypos_new_name)) 
      {
        # print('** FROM')
        old.name = hypos_new_name[from]
      
        # print(old.name)
        idx.from = which(rownames(consensus) == old.name)
        # print(rownames(consensus)[idx.from])
        from = idx.from

      }

   if(to %in% names(hypos_new_name)) 
      {
      # print('** TO')
        old.name = hypos_new_name[to]
      
        # print(old.name)
        idx.to = which(colnames(consensus) == old.name)
        # print(rownames(consensus)[idx.from])
        to = idx.to
      }
        # print(names(hypos_new_name))
      
      ### BUG HERE
     # print('**')
     # print(e)
     
        if(!all(c(from,to) %in% colnames(consensus))) 
        {   
            eAttrs$lwd[e] = MIN.HITS
            eAttrs$label[e] = paste0('', MIN.HITS)
        }
            else 
            {
                    eAttrs$lwd[e] = consensus[from, to]
                eAttrs$label[e] = paste0('', consensus[from, to])
        }
     # sprint(colnames(consensus))
     # print('--')
     # print((hypos_new_name))
     # eAttrs$lwd[e] = consensus[from, to]
     
    }
    # print(sort(colnames(consensus)))

  # remove arrows from logic node (hidden and)
for(e in edge_names) {
    edge = unlist(strsplit(e, '~'))
    from = edge[1]
    to = edge[2]
    
    # print(e)
    # print(is.logic.node.down(to))
    # print(is.logic.node.up(from))
        
    if (is.logic.node.down(to)) {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
     
    # print(from)
    # print(substr(to, start=1, stop=2))
     # eAttrs$color[e] = 'red'
     
     if(substr(to, start=1, stop=2) == 'OR')
         eAttrs$color[e] = 'orange'
     if(substr(to, start=1, stop=3) == 'XOR')
         eAttrs$color[e] = 'red'
     if(substr(to, start=1, stop=3) == 'AND')
         eAttrs$color[e] = 'darkgreen'
           
     eAttrs$lty[e] = 'dashed'

    nAttrs$shape[to] = 'circle' 
    # nAttrs$shape[to] = function(w, ...) {
        # for(i in 1:3) points(w, cex=i, ...)
    # }
    
    # print(nAttrs)

    }

    if (is.logic.node.up(from)) {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
      # eAttrs$color[e] = 'black'
      
      eAttrs$lty[e] = 'dashed'
      
        # print(substr(to, start=1, stop=2))
        # print(substr(from, start=1, stop=2))
    # print(from)
    # print(to)
 
 
      if(substr(from, start=1, stop=4) == 'UPOR')
         eAttrs$color[e] = 'orange'
     if(substr(from, start=1, stop=5) == 'UPXOR')
         eAttrs$color[e] = 'red'
     if(substr(from, start=1, stop=5) == 'UPAND')
         eAttrs$color[e] = 'darkgreen'
     
      
    } else if(substr(from, start=1, stop=1) == '*') {
      eAttrs$logic[e] = T
      eAttrs$arrowsize[e] = 0
      eAttrs$color[e] = 'black'
    }
  }
  
    plot(graph, nodeAttrs=nAttrs, edgeAttrs=eAttrs, main=title, ... )
  
  # Adds the legend to the plot
  # if (legend) {
    # valid_events = colnames(hypo_mat)[which(colnames(hypo_mat) %in% colnames(adjacency))]
    # legend_names = unique(all.events[which(rownames(all.events) %in% valid_events), 'type'])
    # pt_bg = all.types[legend_names, 'color']
    # legend_colors = rep('black', length(legend_names))
    # pch = rep(21, length(legend_names))
    
    # if (length(legend_logic) > 0) {
      # pch = c(pch, 0, 0, rep(22, length(legend_logic)))
      # legend_names = c(legend_names, ' ', expression(bold('Patterns')), names(legend_logic))
      # legend_colors = c(legend_colors, 'white', 'white', rep('black', length(legend_logic)))
      # pt_bg = c(pt_bg, 'white', 'white', legend_logic)  
    # }

    # if (length(legend_pathways) > 0) {
      # pch = c(pch, 0, 0, rep(21, length(legend_pathways)))
      # legend_names = c(legend_names, ' ', expression(bold('Pathways')), names(legend_pathways))
      # pt_bg = c(pt_bg, 'white', 'white', rep('white', length(legend_pathways)))
      # legend_colors = c(legend_colors, 'white', 'white', legend_pathways)  
    # }
    
    # if (legend.pos == 'bottom') {
      # legend.pos.l = 'bottomleft'
      # legend.pos.r = 'bottomright'
    # } else if (legend.pos == 'top') {
      # legend.pos.l = 'topleft'
      # legend.pos.r = 'topright'
    # } else {
      # legend.pos.l = locator(1)
      # legend.pos.r = locator(1)
    # }

    # legend(legend.pos.r,
           # legend = legend_names,
           # title = expression(bold('Events type')),
           # bty = 'n',
           # cex = legend.cex,
           # pt.cex = 1.5 * legend.cex,
           # pch = pch,
           # col = legend_colors,
           # pt.bg = pt_bg)

    # #add thickness legend
    # valid_names = node_names
    # if(expand) {
      # valid_names = node_names[unlist(lapply(node_names, function(x){!is.logic.node(x)}))]
    # }
    # valid_names = grep('^[*]_(.+)$', valid_names, value = T, invert=T)
    # dim = nAttrs$height[valid_names]
    # prob = marginal_p[valid_names, ]
    
    # min = min(dim)
    # p_min = round(min(prob) * 100, 0)
    # max = max(dim)
    # p_max = round(max(prob) * 100, 0)
    
    
    
     # # if( 'Hypothesis' %in% all.types ) hypo.names = rownames(all.events$types='Hypothesis'))
  # # else hypo.names = NA
  
    # nonhypo.names = setdiff(rownames(all.events, hypo.names)
    
    # marginal_p = marginal_p[nonhypo.names, , drop = FALSE]

    # # Get label of the (first) event with minimum marginale 
    # min.p =   rownames(marginal_p)[which(min(marginal_p) == marginal_p) ]
    # label.min = as.events(x)[ min.p[1] , , drop = FALSE]

    # # Get label of the (first) event with max marginale 
    # max.p = rownames(marginal_p)[which(max(marginal_p) == marginal_p) ]
    # label.max = as.events(x)[ max.p[1] , ,  drop = FALSE]

    # # Frequency labels
    # min.freq = round(min(marginal_p) * 100, 0)
    # max.freq = round(max(marginal_p) * 100, 0)
    
    # freq.labels = c( 
      # paste0(min.freq, ifelse((min.freq < 10 && max.freq > 9), '%  ', '%'), ' ', label.min[, 'event'], ' (min)'),
      # paste0(max.freq, '% ', label.max[, 'event'], ' (max)')
    # )
  
    # stat.pch = c(21, 21)
    # pt.bg = c(
        # as.colors(x)[label.min[, 'type']], 
      # as.colors(x)[label.max[, 'type']]
      # )
    # col = c('black', 'black')
        
    # # Further stats
    # y = x
    # if('Hypothesis' %in% as.types(x)) 
            # y = delete.type(x, 'Hypothesis')
            
    # freq.labels = c(freq.labels, 
      # ' ',
      # expression(bold('Sample size')),
      # paste0('n = ', nsamples(y)),
      # paste0('m = ', nevents(y)),
      # paste0('|G| = ', ngenes(y))     
    # ) 
     
    # stat.pch = c(stat.pch, 0, 0, 20,  20, 20)
    # pt.bg = c(pt.bg, 'white', 'white', rep('black', 3))
    # col = c(col, 'white', 'white', rep('black', 3)) 
    
    # legend(legend.pos.l,
           # legend = freq.labels,
           # title = expression(bold('Events frequency')),
           # bty = 'n',
           # box.lty = 3,
           # box.lwd = .3,
           # pch = stat.pch,
           # pt.cex = 1.5  * legend.cex,
           # ncol = 1,
           # pt.bg = pt.bg,
           # cex = legend.cex,
           # col = col)
       
  # }
  
  if(!is.na(file))
  {
    dev.copy2pdf(file = file)
  }
  # cat('\n')
}


