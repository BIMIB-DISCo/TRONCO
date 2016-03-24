#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution


#' Reconstruct a progression model using CAPRESE algorithm
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.caprese(test_dataset_no_hypos)
#'
#' @title tronco caprese
#' @param data A TRONCO compliant dataset.
#' @param lambda Coefficient to combine the raw estimate with a correction factor into a shrinkage estimator. 
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.caprese
#' @importFrom stats phyper
#' @importFrom bnlearn empty.graph set.arc
#' 
tronco.caprese <- function(data,
                           lambda = 0.5,
                           silent = FALSE ) {

    ## Check for the inputs to be correct.
    
    if (is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    if (lambda < 0 || lambda > 1) {
        stop("The value of the shrinkage parameter lambda has to be in [0:1]!",
             call. = FALSE);
    }

    ## check for the input to be compliant.
    
    is.compliant(data)

    ## check if there are hypotheses

    if (npatterns(data) > 0) {
        warning("Patters found in input for tronco.caprese\n")
    }

    ## Reconstruct the reconstruction with CAPRESE.
    
    if (silent == FALSE) {
        cat('*** Checking input events.\n')
        invalid = consolidate.data(data, TRUE)      
        if (length(unlist(invalid)) > 0)
            warning("Input events should be consolidated - see consolidate.data.");

        cat(paste0(
            '*** Inferring a progression model with the following settings.\n',
            '\tDataset size: n = ', nsamples(data), ', m = ', nevents(data), '.\n',
            '\tAlgorithm: CAPRESE with shrinkage coefficient: ', lambda, '.\n'
        ))
    }
    reconstruction = caprese.fit(dataset = data$genotypes,
                                 lambda = lambda,
                                 silent = silent);

    rownames(reconstruction$confidence) =
        c("temporal priority",
          "probability raising",
          "hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";

        ## Set rownames and colnames to the parents positions.
        
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";

        ## Set rownames and colnames to the adjacency matrices.
        
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);

    }

    ## Structure to save the results.
    
    results = data;
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;

    if (!silent) {
        cat('*** Evaluating LogLik informations.\n')

    }

    bayes.net = as.bnlearn.network(results, model = 'caprese')
    logLik = logLik(bayes.net$net, data = bayes.net$data)
    results$model$caprese$logLik = logLik

    ## The reconstruction has been completed.
    if (!silent)
        cat(paste(
            "The reconstruction has been successfully completed in", 
            format(.POSIXct(round(reconstruction$execution.time[3],
                                  digits = 0),
                            tz = "GMT"),
                   "%Hh:%Mm:%Ss"), 
            "\n"));
    
    return(results);
}


#' Reconstruct a progression model using CAPRI algorithm
#'
#' @examples
#' data(test_dataset)
#' recon = tronco.capri(test_dataset, nboot = 5)
#'
#' @title tronco capri
#' @param data A TRONCO compliant dataset.
#' @param command Parameter to define to heuristic search to be performed. Hill Climbing and Tabu search are currently available.
#' @param regularization Select the regularization for the likelihood estimation, e.g., BIC, AIC. 
#' @param do.boot A parameter to disable/enable the estimation of the error rates give the reconstructed model.
#' @param nboot Number of bootstrap sampling (with rejection) to be performed when estimating the selective advantage scores. 
#' @param pvalue Pvalue to accept/reject the valid selective advantage relations. 
#' @param min.boot Minimum number of bootstrap sampling to be performed. 
#' @param min.stat A parameter to disable/enable the minimum number of bootstrap sampling required besides nboot if any sampling is rejected. 
#' @param boot.seed Initial seed for the bootstrap random sampling.
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.capri
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths
#' @importFrom stats phyper AIC BIC wilcox.test
#' 
tronco.capri <- function(data,
                         command = "hc",
                         regularization = c("bic","aic"), 
                         do.boot = TRUE, 
                         nboot = 100, 
                         pvalue = 0.05, 
                         min.boot = 3, 
                         min.stat = TRUE, 
                         boot.seed = NULL, 
                         silent = FALSE ) {

    ## Enforce data to be numeric
    data = enforce.numeric(data)

    ## Check for the inputs to be correct.
    
    if (is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    
    if (is.null(data$hypotheses)) {
        data$hypotheses = NA;
    }
    
    if (command != "hc" && command != "tabu") {
        stop("The inference can be performed either by hill climbing or tabu search!",call. = FALSE);
    }

    if (pvalue < 0 || pvalue > 1) {
        stop("The value of the pvalue has to be in [0:1]!",call. = FALSE);
    }

    if (! all(regularization %in% c('loglik', 'bic', 'aic'))) {
        stop("Possible regularization are loglik, bic or aic",call. = FALSE);
    }

    ## Check for the input to be compliant.
    
    is.compliant(data)

    ## Reconstruct the reconstruction with CAPRI.
    
    if (is.null(boot.seed)) {
        my.seed = "NULL"    
    }
    else {
        my.seed = boot.seed;
    }
    if (silent == FALSE) {
        cat('*** Checking input events.\n')
        invalid = consolidate.data(data, TRUE)
        if (length(unlist(invalid)) > 0) warning(
            "Input events should be consolidated - see consolidate.data."
        );

        cat(paste0(
            '*** Inferring a progression model with the following settings.\n',
            '\tDataset size: n = ',
            nsamples(data),
            ', m = ',
            nevents(data), '.\n',
            '\tAlgorithm: CAPRI with \"',
            paste0(regularization, collapse = ", "),
            '\" regularization and \"',
            command,
            '\" likelihood-fit strategy.\n',
            '\tRandom seed: ',
            my.seed, '.\n',
            '\tBootstrap iterations (Wilcoxon): ',
            ifelse(do.boot, nboot, 'disabled'), '.\n',
            ifelse(do.boot, 
                   paste0('\t\texhaustive bootstrap: ',
                          min.stat,
                          '.\n\t\tp-value: ',
                          pvalue,
                          '.\n\t\tminimum bootstrapped scores: ',
                          min.boot, '.\n'), '')        
        ))
    }

    reconstruction =
        capri.fit(data$genotypes,
                  data$hypotheses,
                  command = command,
                  regularization = regularization,
                  do.boot = do.boot,
                  nboot = nboot,
                  pvalue = pvalue,
                  min.boot = min.boot,
                  min.stat = min.stat,
                  boot.seed = boot.seed,
                  silent = silent);

    rownames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);
    colnames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);

    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";

        ## Set rownames and colnames to the parents positions.
        
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";

        ## Set rownames and colnames to the adjacency matrices.
        
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);

    }

    ## Structure to save the results.
    
    results = data;
    results$adj.matrix.prima.facie = reconstruction$adj.matrix.prima.facie
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;

    ## Add BIC/AIC/LogLik informations

    if (!silent) {
        cat('*** Evaluating BIC / AIC / LogLik informations.\n')
    }

    if ("loglik" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'capri_loglik')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$capri_loglik$score = score
        results$model$capri_loglik$logLik = logLik
    }

    if ("bic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'capri_bic')
        score = BIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$capri_bic$score = score
        results$model$capri_bic$logLik = logLik
    }

    if ("aic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'capri_aic')
        score = AIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$capri_aic$score = score
        results$model$capri_aic$logLik = logLik
    }

    ## the reconstruction has been completed.
    
    if (!silent)
        cat(paste(
            "The reconstruction has been successfully completed in", 
            format(.POSIXct(round(reconstruction$execution.time[3],
                                  digits = 0),
                            tz = "GMT"),
                   "%Hh:%Mm:%Ss"), 
            "\n"));

    return(results);
}


#' Reconstruct a progression model using Edmonds algorithm combined 
#' with probabilistic causation
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.mst.edmonds(test_dataset_no_hypos, nboot = 5)
#'
#' @title tronco mst edmonds
#' @param data A TRONCO compliant dataset.
#' @param regularization Select the regularization for the 
#' likelihood estimation, e.g., BIC, AIC. 
#' @param do.boot A parameter to disable/enable the estimation 
#' of the error rates give the reconstructed model.
#' @param nboot Number of bootstrap sampling (with rejection) 
#' to be performed when estimating the selective advantage scores. 
#' @param pvalue Pvalue to accept/reject the valid selective 
#' advantage relations. 
#' @param min.boot Minimum number of bootstrap sampling to be 
#' performed. 
#' @param min.stat A parameter to disable/enable the minimum number
#' of bootstrap sampling required besides nboot if any sampling 
#' is rejected. 
#' @param boot.seed Initial seed for the bootstrap random sampling.
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.mst.edmonds
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths
#' @importFrom infotheo mutinformation
#' @importFrom stats phyper AIC BIC
#' 
tronco.mst.edmonds <- function(data,
                               regularization = "no_reg", 
                               do.boot = TRUE, 
                               nboot = 100, 
                               pvalue = 0.05, 
                               min.boot = 3, 
                               min.stat = TRUE, 
                               boot.seed = NULL, 
                               silent = FALSE ) {

    ## Enforce data to be numeric
    data = enforce.numeric(data)

    ## Check for the inputs to be correct.
    
    if (is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    
    if (is.null(data$hypotheses)) {
        data$hypotheses = NA;
    }
    
    if (pvalue < 0 || pvalue > 1) {
        stop("The value of the pvalue has to be in [0:1]!",call. = FALSE);
    }

    if (! all(regularization %in% c('no_reg', 'loglik', 'bic', 'aic'))) {
        stop("Possible regularization are no-reg, loglik, bic or aic",call. = FALSE);
    }

    ## Check for the input to be compliant.
    
    is.compliant(data)

    ## check if there are hypotheses

    if (npatterns(data) > 0) {
        warning("Patters found in input for tronco.mst.edmonds\n")
    }

    ## Reconstruct the reconstruction with Edmonds.
    
    if (is.null(boot.seed)) {
        my.seed = "NULL"    
    }
    else {
        my.seed = boot.seed;
    }
    if (silent == FALSE) {
        cat('*** Checking input events.\n')
        invalid = consolidate.data(data, TRUE)
        if (length(unlist(invalid)) > 0) warning(
            "Input events should be consolidated - see consolidate.data."
        );

        cat(paste0(
            '*** Inferring a progression model with the following settings.\n',
            '\tDataset size: n = ',
            nsamples(data),
            ', m = ',
            nevents(data), '.\n',
            '\tAlgorithm: Edmonds with \"',
            paste0(regularization, collapse = ", "),
            '\" regularization',
            '\tRandom seed: ',
            my.seed, '.\n',
            '\tBootstrap iterations (Wilcoxon): ',
            ifelse(do.boot, nboot, 'disabled'), '.\n',
            ifelse(do.boot, 
                   paste0('\t\texhaustive bootstrap: ',
                          min.stat,
                          '.\n\t\tp-value: ',
                          pvalue,
                          '.\n\t\tminimum bootstrapped scores: ',
                          min.boot, '.\n'), '')        
        ))
    }

    reconstruction =
        edmonds.fit(data$genotypes,
                    regularization = regularization,
                    do.boot = do.boot,
                    nboot = nboot,
                    pvalue = pvalue,
                    min.boot = min.boot,
                    min.stat = min.stat,
                    boot.seed = boot.seed,
                    silent = silent)

    rownames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes)
    colnames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes)

    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test")
    colnames(reconstruction$confidence) = "confidence"
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes)
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes)
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes)
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes)
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes)
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes)

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability"
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes)
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability"

        ## Set rownames and colnames to the parents positions.
        
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$parents.pos) = "parents"

        ## Set rownames and colnames to the adjacency matrices.
        
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes)
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes)
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes)

    }

    ## Structure to save the results.
    
    results = data
    results$adj.matrix.prima.facie = reconstruction$adj.matrix.prima.facie
    results$confidence = reconstruction$confidence
    results$model = reconstruction$model
    results$parameters = reconstruction$parameters
    results$execution.time = reconstruction$execution.time

    ## Add BIC/AIC/LogLik informations

    if (!silent) {
        cat('*** Evaluating BIC / AIC / LogLik informations.\n')
    }

    if ("no_reg" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'edmonds_no_reg')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$edmonds_no_reg$score = score
        results$model$edmonds_no_reg$logLik = logLik
    }
    
    if ("loglik" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'edmonds_loglik')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$edmonds_loglik$score = score
        results$model$edmonds_loglik$logLik = logLik
    }

    if ("bic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'edmonds_bic')
        score = BIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$edmonds_bic$score = score
        results$model$edmonds_bic$logLik = logLik
    }

    if ("aic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'edmonds_aic')
        score = AIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$edmonds_aic$score = score
        results$model$edmonds_aic$logLik = logLik
    }

    ## the reconstruction has been completed.
    
    if (!silent)
        cat(paste(
            "The reconstruction has been successfully completed in", 
            format(.POSIXct(round(reconstruction$execution.time[3],
                                  digits = 0),
                            tz = "GMT"),
                   "%Hh:%Mm:%Ss"), 
            "\n"));

    return(results);
}


#' Reconstruct a progression model using Chow Liu 
#' algorithm combined with probabilistic causation
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.mst.chowliu(test_dataset_no_hypos, nboot = 5)
#'
#' @title tronco mst chow liu
#' @param data A TRONCO compliant dataset.
#' @param regularization Select the regularization for the
#' likelihood estimation, e.g., BIC, AIC. 
#' @param do.boot A parameter to disable/enable the estimation 
#' of the error rates give the reconstructed model.
#' @param nboot Number of bootstrap sampling (with rejection) 
#' to be performed when estimating the selective advantage scores. 
#' @param pvalue Pvalue to accept/reject the valid selective
#' advantage relations. 
#' @param min.boot Minimum number of bootstrap sampling 
#' to be performed. 
#' @param min.stat A parameter to disable/enable the minimum 
#' number of bootstrap sampling required besides nboot if 
#' any sampling is rejected. 
#' @param boot.seed Initial seed for the bootstrap random 
#' sampling.
#' @param silent A parameter to disable/enable verbose 
#' messages.
#' @return A TRONCO compliant object with reconstructed 
#' model
#' @export tronco.mst.chowliu
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths
#' @importFrom gRapHD minForest
#' @importFrom stats phyper AIC BIC
#' 
tronco.mst.chowliu <- function(data,
                               regularization = c("bic","aic"), 
                               do.boot = TRUE, 
                               nboot = 100, 
                               pvalue = 0.05, 
                               min.boot = 3, 
                               min.stat = TRUE, 
                               boot.seed = NULL, 
                               silent = FALSE ) {

    ## Enforce data to be numeric
    data = enforce.numeric(data)

    ## Check for the inputs to be correct.
    
    if (is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    
    if (is.null(data$hypotheses)) {
        data$hypotheses = NA;
    }
    
    if (pvalue < 0 || pvalue > 1) {
        stop("The value of the pvalue has to be in [0:1]!",call. = FALSE);
    }

    if (! all(regularization %in% c('loglik','bic', 'aic'))) {
        stop("Possible regularization are loglik, bic or aic",call. = FALSE);
    }

    ## Check for the input to be compliant.
    
    is.compliant(data)

    ## check if there are hypotheses

    if (npatterns(data) > 0) {
        warning("Patters found in input for tronco.chow.liu\n")
    }

    ## Reconstruct the reconstruction with Chow Liu.
    
    if (is.null(boot.seed)) {
        my.seed = "NULL"    
    }
    else {
        my.seed = boot.seed;
    }
    if (silent == FALSE) {
        cat('*** Checking input events.\n')
        invalid = consolidate.data(data, TRUE)
        if (length(unlist(invalid)) > 0) warning(
            "Input events should be consolidated - see consolidate.data."
        );

        cat(paste0(
            '*** Inferring a progression model with the following settings.\n',
            '\tDataset size: n = ',
            nsamples(data),
            ', m = ',
            nevents(data), '.\n',
            '\tAlgorithm: Chow Liu with \"',
            paste0(regularization, collapse = ", "),
            '\" regularization',
            '\tRandom seed: ',
            my.seed, '.\n',
            '\tBootstrap iterations (Wilcoxon): ',
            ifelse(do.boot, nboot, 'disabled'), '.\n',
            ifelse(do.boot, 
                   paste0('\t\texhaustive bootstrap: ',
                          min.stat,
                          '.\n\t\tp-value: ',
                          pvalue,
                          '.\n\t\tminimum bootstrapped scores: ',
                          min.boot, '.\n'), '')        
        ))
    }

    reconstruction =
        chow.liu.fit(data$genotypes,
                  regularization = regularization,
                  do.boot = do.boot,
                  nboot = nboot,
                  pvalue = pvalue,
                  min.boot = min.boot,
                  min.stat = min.stat,
                  boot.seed = boot.seed,
                  silent = silent);

    rownames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);
    colnames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);

    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";

        ## Set rownames and colnames to the parents positions.
        
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";

        ## Set rownames and colnames to the adjacency matrices.
        
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);

    }

    ## Structure to save the results.
    
    results = data;
    results$adj.matrix.prima.facie = reconstruction$adj.matrix.prima.facie
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;

    ## Add BIC/AIC/LogLik informations

    if (!silent) {
        cat('*** Evaluating BIC / AIC / LogLik informations.\n')
    }

    if ("loglik" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'chow_liu_loglik')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$chow_liu_loglik$score = score
        results$model$chow_liu_loglik$logLik = logLik
    }

    if ("bic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'chow_liu_bic')
        score = BIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$chow_liu_bic$score = score
        results$model$chow_liu_bic$logLik = logLik
    }

    if ("aic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'chow_liu_aic')
        score = AIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$chow_liu_aic$score = score
        results$model$chow_liu_aic$logLik = logLik
    }

    ## the reconstruction has been completed.
    
    if (!silent)
        cat(paste(
            "The reconstruction has been successfully completed in", 
            format(.POSIXct(round(reconstruction$execution.time[3],
                                  digits = 0),
                            tz = "GMT"),
                   "%Hh:%Mm:%Ss"), 
            "\n"));

    return(results);
}


#' Reconstruct a progression model using Prim algorithm combined with probabilistic causation
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.mst.prim(test_dataset_no_hypos, nboot = 5)
#'
#' @title tronco mst prim
#' @param data A TRONCO compliant dataset.
#' @param regularization Select the regularization for the 
#' likelihood estimation, e.g., BIC, AIC. 
#' @param do.boot A parameter to disable/enable the estimation 
#' of the error rates give the reconstructed model.
#' @param nboot Number of bootstrap sampling (with rejection) 
#' to be performed when estimating the selective advantage scores. 
#' @param pvalue Pvalue to accept/reject the valid selective 
#' advantage relations. 
#' @param min.boot Minimum number of bootstrap sampling to 
#' be performed. 
#' @param min.stat A parameter to disable/enable the minimum 
#' number of bootstrap sampling required besides nboot if 
#' any sampling is rejected. 
#' @param boot.seed Initial seed for the bootstrap random sampling.
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.mst.prim
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph get.edgelist E E<-
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths minimum.spanning.tree
#' @importFrom infotheo mutinformation
#' @importFrom stats phyper AIC BIC
#' 
tronco.mst.prim <- function(data,
                            regularization = "no_reg", 
                            do.boot = TRUE, 
                            nboot = 100, 
                            pvalue = 0.05, 
                            min.boot = 3, 
                            min.stat = TRUE, 
                            boot.seed = NULL, 
                            silent = FALSE ) {

    ## Enforce data to be numeric
    data = enforce.numeric(data)

    ## Check for the inputs to be correct.
    
    if (is.null(data) || is.null(data$genotypes)) {
        stop("The dataset given as input is not valid.");
    }
    
    if (is.null(data$hypotheses)) {
        data$hypotheses = NA;
    }
    
    if (pvalue < 0 || pvalue > 1) {
        stop("The value of the pvalue has to be in [0:1]!",call. = FALSE);
    }

    if (! all(regularization %in% c('no_reg','loglik', 'bic', 'aic'))) {
        stop("Possible regularization are no-reg, loglik, bic or aic",call. = FALSE);
    }

    ## Check for the input to be compliant.
    
    is.compliant(data)

    ## check if there are hypotheses

    if (npatterns(data) > 0) {
        warning("Patters found in input for tronco.prim\n")
    }

    ## Reconstruct the reconstruction with Prim.
    
    if (is.null(boot.seed)) {
        my.seed = "NULL"    
    }
    else {
        my.seed = boot.seed;
    }
    if (silent == FALSE) {
        cat('*** Checking input events.\n')
        invalid = consolidate.data(data, TRUE)
        if (length(unlist(invalid)) > 0) warning(
            "Input events should be consolidated - see consolidate.data."
        );

        cat(paste0(
            '*** Inferring a progression model with the following settings.\n',
            '\tDataset size: n = ',
            nsamples(data),
            ', m = ',
            nevents(data), '.\n',
            '\tAlgorithm: Prim with \"',
            paste0(regularization, collapse = ", "),
            '\" regularization',
            '\tRandom seed: ',
            my.seed, '.\n',
            '\tBootstrap iterations (Wilcoxon): ',
            ifelse(do.boot, nboot, 'disabled'), '.\n',
            ifelse(do.boot, 
                   paste0('\t\texhaustive bootstrap: ',
                          min.stat,
                          '.\n\t\tp-value: ',
                          pvalue,
                          '.\n\t\tminimum bootstrapped scores: ',
                          min.boot, '.\n'), '')        
        ))
    }

    reconstruction =
        prim.fit(data$genotypes,
                 regularization = regularization,
                 do.boot = do.boot,
                 nboot = nboot,
                 pvalue = pvalue,
                 min.boot = min.boot,
                 min.stat = min.stat,
                 boot.seed = boot.seed,
                 silent = silent);

    rownames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);
    colnames(reconstruction$adj.matrix.prima.facie) = colnames(data$genotypes);

    rownames(reconstruction$confidence) = c("temporal priority","probability raising","hypergeometric test");
    colnames(reconstruction$confidence) = "confidence";
    rownames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[1,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[2,1]]) = colnames(data$genotypes);
    rownames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);
    colnames(reconstruction$confidence[[3,1]]) = colnames(data$genotypes);

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability";
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability";

        ## Set rownames and colnames to the parents positions.
        
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$parents.pos) = "parents";

        ## Set rownames and colnames to the adjacency matrices.
        
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.pf) = colnames(data$genotypes);
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(data$genotypes);

    }

    ## Structure to save the results.
    
    results = data;
    results$adj.matrix.prima.facie = reconstruction$adj.matrix.prima.facie
    results$confidence = reconstruction$confidence;
    results$model = reconstruction$model;
    results$parameters = reconstruction$parameters;
    results$execution.time = reconstruction$execution.time;

    ## Add BIC/AIC/LogLik informations

    if (!silent) {
        cat('*** Evaluating BIC / AIC / LogLik informations.\n')
    }

    if ("no_reg" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'prim_no_reg')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$prim_no_reg$score = score
        results$model$prim_no_reg$logLik = logLik
    }

    if ("loglik" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'prim_no_loglik')
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
        results$model$prim_loglik$score = score
        results$model$prim_loglik$logLik = logLik
    }

    if ("bic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'prim_bic')
        score = BIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$prim_bic$score = score
        results$model$prim_bic$logLik = logLik
    }

    if ("aic" %in% regularization) {
        bayes.net = as.bnlearn.network(results, model = 'prim_aic')
        score = AIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
        results$model$prim_aic$score = score
        results$model$prim_aic$logLik = logLik
    }

    ## the reconstruction has been completed.
    
    if (!silent)
        cat(paste(
            "The reconstruction has been successfully completed in", 
            format(.POSIXct(round(reconstruction$execution.time[3],
                                  digits = 0),
                            tz = "GMT"),
                   "%Hh:%Mm:%Ss"), 
            "\n"));

    return(results);
}


#' Bootstrap a reconstructed progression model
#'
#' @examples
#' data(test_model)
#' boot = tronco.bootstrap(test_model, nboot = 5)
#'
#' @title tronco bootstrap
#' @param reconstruction The output of tronco.capri or 
#' tronco.caprese
#' @param type Parameter to define the type of sampling 
#' to be performed, e.g., non-parametric for uniform sampling.
#' @param nboot Number of bootstrap sampling to be performed 
#' when estimating the model confidence.
#' @param verbose Should I be verbose?
#' @param cores.ratio Percentage of cores to use
#' coresRate * (numCores - 1)
#' @return A TRONCO compliant object with reconstructed model
#' @importFrom doParallel registerDoParallel  
#' @importFrom foreach foreach %dopar%
#' @importFrom iterators icount
#' @importFrom parallel stopCluster makeCluster detectCores
#' @export tronco.bootstrap
#' 
tronco.bootstrap <- function(reconstruction,
                             type = "non-parametric",
                             nboot = 100,
                             verbose = FALSE,
                             cores.ratio = 1) {
    
    ## Check for the inputs to be given.
    
    if (is.null(reconstruction)) {
        stop("A valid reconstruction has to be provided in order to estimate its confidence.",
             call. = FALSE)
    }

    ## Check for the input to be compliant.
    
    is.compliant(reconstruction)
    
    if (type == "statistical"
        && !((reconstruction$parameters$algorithm == "CAPRI"
             || reconstruction$parameters$algorithm == "PRIM"
             || reconstruction$parameters$algorithm == "CHOW_LIU"
             || reconstruction$parameters$algorithm == "EDMONDS")
             && reconstruction$parameters$do.boot == TRUE)) {
        stop(paste("To perform statistical bootstrap, the algorithm used for",
                   "the reconstruction must be CAPRI, PRIM, CHOW_LIU or EDMONDS",
                   "with bootstrap."),
             call. = FALSE)
    }

    ## Set all the needed parameters to perform the bootstrap
    ## estimation.

    if (!type %in% c("non-parametric", "statistical")) {
        stop(paste("The types of bootstrap that can be performed are:",
                   "non-parametric or statistical."),
             call. = FALSE)
    }
    
    ## Perform the selected bootstrap procedure.
    
    cat("*** Executing now the bootstrap procedure, this may take a long time...\n")

    parameters = as.parameters(reconstruction)

    if (parameters$algorithm == "CAPRESE") {

        lambda = parameters$lambda
        curr.boot = bootstrap(reconstruction, 
                              type,
                              nboot,
                              cores.ratio)
        cat("Performed", type,
            "bootstrap with", nboot,
            "resampling and", lambda, 
            "as shrinkage parameter.\n")

    } else {

        curr.boot = bootstrap(reconstruction, 
                              type,
                              nboot,
                              cores.ratio)

        cat("Performed", type,
            "bootstrap with", nboot,
            "resampling")

        if (parameters$do.boot == TRUE) {
            cat(" and", 
                parameters$pvalue,
                "as pvalue for the statistical tests")
        } 
        cat(".\n")
    }
    reconstruction$bootstrap = curr.boot
    return(reconstruction)
}


#' Plots a progression model from a recostructed dataset
#' @title tronco.plot
#'
#' @examples
#' data(test_model)
#' tronco.plot(test_model)
#'
#' @param x A reconstructed model (the output of the inference 
#' by a tronco function)
#' @param models A vector containing the names of the 
#' algorithms used (caprese, capri_bic, etc)
#' @param fontsize For node names. Default NA for 
#' automatic rescaling
#' @param height Proportion node height - node width. 
#' Default height 2
#' @param width Proportion node height - node width. 
#' Default width 2
#' @param height.logic Height of logical nodes. 
#' Defaul 1
#' @param pf Should I print Prima Facie? 
#' Default False
#' @param disconnected Should I print disconnected 
#' nodes? Default False
#' @param scale.nodes Node scaling coefficient 
#' (based on node frequency). Default NA (autoscale)
#' @param title Title of the plot. Default as.description(x)  
#' @param confidence Should I add confidence 
#' informations? No if NA
#' @param p.min p-value cutoff. Default automatic
#' @param legend Should I visualise the legend?
#' @param legend.cex CEX value for legend. Default 1.0
#' @param edge.cex CEX value for edge labels. Default 1.0
#' @param label.edge.size Size of edge labels. 
#' Default NA for automatic rescaling
#' @param expand Should I expand hypotheses? Default TRUE
#' @param genes Visualise only genes in this list. 
#' Default NULL, visualise all.
#' @param relations.filter Filter relations to dispaly 
#' according to this functions. Default NA
#' @param edge.color Edge color. Default 'black'
#' @param pathways.color RColorBrewer colorser 
#' for patways. Default 'Set1'.
#' @param file String containing filename for PDF output. 
#' If NA no PDF output will be provided
#' @param legend.pos Legend position. Default 'bottom',
#' @param pathways A vector containing pathways information 
#' as described in as.patterns()
#' @param lwd Edge base lwd. Default 3
#' @param samples.annotation = List of samples to search 
#' for events in model
#' @param export.igraph If TRUE export the igraph 
#' object generated
#' @param ... Additional arguments for RGraphviz 
#' plot function
#' @return Information about the reconstructed model               
#' @export tronco.plot
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph V V<- igraph.to.graphNEL igraph.from.graphNEL
#' @importFrom Rgraphviz edgeNames plot
#' @importFrom graphics locator
#' @importFrom grDevices dev.copy2pdf
#' 
tronco.plot <- function(x,
                        models = names(x$model),
                        fontsize = NA, 
                        height=2,
                        width=3,
                        height.logic = 1,
                        pf = FALSE, 
                        disconnected=FALSE,
                        scale.nodes=NA,
                        title = as.description(x),  
                        confidence = NA, 
                        p.min = 0.05,
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
                        samples.annotation = NA,
                        export.igraph = FALSE,
                        ...) {
    hidden.and = FALSE

    ## Checks if reconstruction exists.
    
    if (missing(x)) {
        stop("reconstruction missing, usage: hypo.plot(reconstruction, ...",
             call. = FALSE);
    }

    logical_op = list("AND", "OR", "NOT", "XOR", "*", "UPAND", "UPOR", "UPXOR")


    if (length(models) > 2) {
        stop("Too many regularizators (max is 2)", call. = FALSE)
    }

    if (! all(models %in% names(x$model))) {
        stop(paste(paste(models, collapse=' '), 
                   "not in reconstructed models. Use: ",
                   paste(names(x$model), collapse=' ')),
        call. = FALSE);
    }

    if (!is.na(samples.annotation) && !is.null(pathways))
        stop('Select either to annotate pathways or a sample.')

    ## Annotate samples.
    
    if (!all(is.na(samples.annotation))) {  
        if (!all(samples.annotation %in% as.samples(x))) {
            stop('Sample(s) to annotate are not in the dataset -- see as.samples.')
        }

        if (npatterns(x) > 0) {
            nopatt.data = delete.model(x)
            nopatt.data = delete.type(nopatt.data, 'Pattern')
        } else {
            nopatt.data = x
        }

        sample.events = Reduce(rbind, as.events.in.sample(nopatt.data, samples.annotation))
        sample.events = unique(sample.events[, 'event'])

        cat('Annotating sample',
            samples.annotation,
            'with color red. Annotated genes:',
            paste(sample.events, collapse = ', '),
            '\n')

        pathways = list(sample.events)
        names(pathways) = paste(samples.annotation, collapse = ', ')
        if (nchar(names(pathways)) > 15)
            names(pathways) = paste0(substr(names(pathways), 1, 15), '...')

        pathways.color = 'red'
    }

    sec = ifelse(length(models) == 2, TRUE, FALSE)

    if (sec && !models[2] %in% names(x$model)) {
        stop(paste(models[2], "not in model"), call.=FALSE);
    }

    ## Models objects.
    
    primary = as.models(x, models = models[1])[[1]]    
    if (sec) 
        secondary = as.models(x, models = models[2])[[1]]

    ## USARE getters adj.matrix.
    
    if (sec && !all(rownames(primary$adj.matrix$adj.matrix.fit)
                    %in% rownames(secondary$adj.matrix$adj.matrix.fit))) {     
        stop("primary and secondary must have the same adj.matrix! See: the function tronco.bootstrap.",
             call. = FALSE)
    }

    ## Get the adjacency matrix - this could have been donw with
    ## getters.
    
    adj.matrix = primary$adj.matrix
    if (sec) adj.matrix = secondary$adj.matrix
    c_matrix = adj.matrix$adj.matrix.fit

    if (is.function(relations.filter)) {
        cat('*** Filtering relations according to function "relations.filter", visualizing:\n')
        adj.matrix = as.adj.matrix(x, models = models)
        sel.relation = as.selective.advantage.relations(x, models = models)

        ## Select only relations which get TRUE by "relations.filter".
        
        sel.relation = lapply(sel.relation,
            function(z) {
                ## Apply can not be used - implicit coercion to char is crap
                ## z[ apply(z, 1, relations.filter), ]

                mask = rep(TRUE, nrow(z))
                for(i in 1:nrow(z))
                    mask[i] = relations.filter(z[i, ]) 
                return(z[mask, , drop = FALSE])                              
            })

        sel.relation = get(models[2], sel.relation)

        c_matrix.names = rownames(c_matrix)
        c_matrix = matrix(0, nrow = nrow(c_matrix), ncol = ncol(c_matrix))
        rownames(c_matrix) = c_matrix.names
        colnames(c_matrix) = c_matrix.names

        cat(paste0('Selected relations: ',
                   nrow(sel.relation),
                   ' [out of ',
                   nrow(as.selective.advantage.relations(x,
                                                         models = models)[[2]]),
                   ']\n'))

        if (nrow(sel.relation) > 0) {
            for(i in 1:nrow(sel.relation)) {
                c_matrix[nameToKey(x, sel.relation[i, 'SELECTS']),
                         nameToKey(x, sel.relation[i, 'SELECTED'])] = 1
            }
        }


    }

    ## Get the probabilities.
    
    probabilities = primary$probabilities
    if (sec) {
        probabilities = secondary$probabilities
    }
    marginal_p = probabilities$probabilities.observed$marginal.probs

    ## If prima facie change the adj matrix.
    if (pf) {
        c_matrix = adj.matrix$adj.matrix.pf
    }

    if (all(c_matrix == FALSE)
        || (sec && all(primary$adj.matrix$adj.matrix.fit == FALSE))) {
        warning('No edge in adjacency matrix! Nothing to show here.')
        return(NULL)
    }

    ## Get algorithm parameters.
    
    parameters = x$parameters

    ## Get hypotheses.
    
    hypotheses = x$hypotheses
    hstruct = NULL
    if (!is.null(hypotheses) && !is.na(hypotheses) ) {
        hstruct = hypotheses$hstructure
    }

    ## Get event from genes list.
    
    events = NULL
    if (is.vector(genes)) {
        events =
            unlist(lapply(genes,
                          function(x) {
                              names(which(as.events(x)[, 'event'] == x))
                          }))
    }

    cat('*** Expanding hypotheses syntax as graph nodes:')

    ## Expand hypotheses.
    
    expansion =
        hypotheses.expansion(c_matrix, 
                             hstruct, 
                             hidden.and, 
                             expand,
                             events)
    hypo_mat = expansion[[1]]
    hypos_new_name = expansion[[2]]

    cat('\n*** Rendering graphics\n')

    ## Remove disconnected nodes.
    
    if (!disconnected) { 
        cat('Nodes with no incoming/outgoing edges will not be displayed.\n')
        del = which(rowSums(hypo_mat) + colSums(hypo_mat) == 0 )
        w = !(rownames(hypo_mat) %in% names(del))
        hypo_mat = hypo_mat[w,]
        hypo_mat = hypo_mat[,w]
    }

    attrs = list(node = list())

    hypo_graph = graph.adjacency(hypo_mat)

    v_names = gsub("_.*$", "", V(hypo_graph)$name)
    if (!expand) {
        v_names = gsub("^[*]_(.+)", "*", V(hypo_graph)$name)
    }
    new_name = list()

    for(v in v_names) {
        if (v %in% rownames(x$annotations)) {
            n = x$annotations[v,"event"]
            new_name = append(new_name, n)
        } else {
            new_name = append(new_name, v)
        }
    }

    V(hypo_graph)$label = new_name

    graph = igraph.to.graphNEL(hypo_graph)

    node_names = V(hypo_graph)$name

    nAttrs = list()

    nAttrs$label = V(hypo_graph)$label
    names(nAttrs$label) = node_names

    ## Set a default color.
    
    nAttrs$fillcolor =  rep('White', length(node_names))
    names(nAttrs$fillcolor) = node_names

    ## Set fontsize.

    if (is.na(fontsize)) {
        fontsize = 24 - 4*log(nrow(hypo_mat))
        cat('Set automatic fontsize scaling for node labels: ', fontsize, '\n')
    }
    nAttrs$fontsize = rep(fontsize, length(node_names))
    names(nAttrs$fontsize) = node_names

    ## Set node shape.
    
    nAttrs$shape = rep('ellipse', length(node_names))
    names(nAttrs$shape) = node_names

    ## Set node height.
    
    nAttrs$height = rep(height, length(node_names))
    names(nAttrs$height) = node_names

    ## Set node width.
    
    nAttrs$width = rep(width, length(node_names))
    names(nAttrs$width) = node_names

    short.label = nAttrs$label
    names(short.label) = names(nAttrs$label)
    if (!is.na(scale.nodes)) {

        ## Foreach node.
        
        min_p = min(marginal_p)
        max_p = max(marginal_p)

        for (node in node_names) {
            prefix = gsub("_.*$", "", node)
            if ( !(prefix %in% logical_op)) {
                
                ## Scaling ANDRE CITROLO.
                
                increase_coeff = scale.nodes + (marginal_p[node,] - min_p) / (max_p - min_p)
                nAttrs$width[node] = nAttrs$width[node] * increase_coeff
                nAttrs$height[node] = nAttrs$height[node] * increase_coeff
                nAttrs$label[node] =
                    paste0(nAttrs$label[node],
                           '\\\n',
                           round(marginal_p[node, ]*100, 0),
                           '%',
                           ' (',
                           sum(as.genotypes(x)[, node]),
                           ')')
            }
        }
    }

    ## Use colors defined in
    ## tronco$types.
    
    w =
        unlist(lapply(names(nAttrs$fillcolor),
                      function(w) {
                          if (w %in% rownames(x$annotations)) {
                              x$types[x$annotations[w, 'type'], 'color']
                          } else
                              'White'
                      }))
    nAttrs$fillcolor[] = w

    legend_logic = NULL

    ## Set color, size form and shape each logic nodes (if hypos
    ## expansion actived)
    
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

    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == '*'
    if (any(w)) {
        legend_logic['Co-occurence'] = 'darkgreen'
    }
    nAttrs$fillcolor[which(w)] = 'darkgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic


    ## Node border to black.
    
    nAttrs$color = rep("black", length(node_names))
    names(nAttrs$color) = node_names

    nAttrs$fontcolor = rep("black", length(node_names))
    names(nAttrs$fontcolor) = node_names

    nAttrs$lwd = rep(1, length(node_names))
    names(nAttrs$lwd) = node_names

    ## Set node border based on pathways information.
    
    legend_pathways = NULL
    
    if (!is.null(pathways)) {
        cat('Annotating nodes with pathway information. \n')

        if (length(pathways.color) == 1
            && pathways.color %in% rownames(brewer.pal.info)) {
            cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
            cols = brewer.pal(n=length(names(pathways)), name=pathways.color)
        } else {
            if (length(pathways.color) != length(names(pathways))) 
                stop('You did not provide enough colors to annotate',
                     length(names(pathways)),
                     'pathways. 
                      Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')
            cols = pathways.color
        }

        names(cols) = names(pathways)
        names(nAttrs$col) = node_names

        for (path in names(pathways)) {
            n = short.label[which(short.label %in% pathways[[path]])]
            nAttrs$color[unlist(names(n))] = cols[[path]]
            nAttrs$fontcolor[unlist(names(n))] = cols[[path]]

            nAttrs$lwd[unlist(names(n))] = 4

            if (length(n) > 0) {
                legend_pathways[path] = cols[[path]]
            }
        }
    }

    ## Edges properties.

    edge_names = edgeNames(graph)
    eAttrs = list()

    ## Set temporary edge shape.
    
    eAttrs$lty = rep("solid", length(edge_names))
    names(eAttrs$lty) = edge_names

    ## Set temporary fontocolor.
    
    eAttrs$fontcolor = rep("darkblue", length(edge_names))
    names(eAttrs$fontcolor) = edge_names

    ## Set edge thikness based on prob.
    
    eAttrs$lwd = rep(1, length(edge_names))
    names(eAttrs$lwd) = edge_names

    ## Set edge name based on prob.
    
    eAttrs$label = rep('', length(edge_names))
    names(eAttrs$label) = edge_names
    
    ## Set fontsize to label.edge.size (default)
    
    if (is.na(label.edge.size)) {
        label.edge.size = fontsize/2      
        cat(paste0('Set automatic fontsize for edge labels: ', label.edge.size, '\n'))    
    }
    eAttrs$fontsize = rep(label.edge.size, length(edge_names))
    names(eAttrs$fontsize) = edge_names

    ## Set edge color to black (default).
    
    eAttrs$color = rep(ifelse(sec, 'darkgrey', edge.color), length(edge_names))
    names(eAttrs$color) = edge_names

    ## Set edge arrowsize to 1 (default).
    
    eAttrs$arrowsize = rep(1 * edge.cex, length(edge_names))
    names(eAttrs$arrowsize) = edge_names

    ## Record logic edge.
    
    eAttrs$logic = rep(FALSE, length(edge_names))
    names(eAttrs$logic) = edge_names

    pval.names = c('hg', 'pr', 'tp')
    boot.names = c('npb', 'pb', 'sb')

    edge_label = function(value, conf, edge, model, pvalue) {
        ret = list()

        if (conf %in% boot.names) {
            ret$lwd = 1 + (value * 3)
        } else {
            ret$lwd = eAttrs$lwd[edge]
        }

        if (c %in% c('posterr', 'prederr')) {
            value = mean(unlist(value))
            ret$label =
                paste0(eAttrs$label[edge],
                       ifelse(value < 0.01,
                              "< 0.01",
                              round(value, 2)))
        } else {
            ret$label =
                paste0(eAttrs$label[edge],
                       ifelse(value < 0.01,
                              "< 0.01",
                              round(value, 2)))
        }


        # insert here edges visualization rules

        if (c == 'hg' && value > pvalue) {
            ret$fontcolor = 'red'
            ret$label = paste0(ret$label, ' *')
        }
        else if (c %in% c('pr', 'tp') 
            && model != 'caprese' 
            && value > pvalue) {
            ret$fontcolor = 'red'
            ret$label = paste0(ret$label, ' *')
        } else {
            ret$fontcolor = eAttrs$fontcolor[edge]
        }
        ret$label = paste0(ret$label, '\\\n')
        return(ret)
    }

    if (any(!is.na(confidence))) {
        cat('Adding confidence information: ')
        conf = as.confidence(x, confidence, models)
        cat(paste(paste(confidence, collapse = ', '), '\n'))


        for (c in confidence) {

            if (c == 'eloss') {
                next
                cat('skip eloss \n')
            }

            if (c == 'posterr') {
                conf_sel = as.kfold.posterr(x, models=models, table = TRUE)
            } else if (c == 'prederr'){
                conf_sel = as.kfold.prederr(x, models=models, table = TRUE)
            } else {
                conf_sel = get(c, as.confidence(x, c))
            }
             
            for(e in edge_names) {

                # configure conf_from and conf_to

                edge = unlist(strsplit(e, '~'))

                from = edge[1]
                to = edge[2]

                if (is.logic.node.up(from) || is.logic.node.down(to)) {
                    next
                }


                if (from %in% names(hypos_new_name)) {
                    conf_from = hypos_new_name[[from]]
                } else {
                    conf_from = from
                }
                if (to %in% names(hypos_new_name)) {
                    conf_to = hypos_new_name[[to]]
                } else {
                    conf_to = to
                }

                conf_p = conf_sel
                mod = ''

                if (!c %in% pval.names) {
                    if (sec && primary$adj.matrix$adj.matrix.fit[conf_from, conf_to] == 0) {
                        conf_p = get(models[[2]], conf_sel)
                        mod = models[[2]]
                    } else {
                        conf_p = get(models[[1]], conf_sel)
                        mod = models[[1]]
                    }
                }

                if (c == 'prederr') {
                    value = conf_p[[conf_to]]
                } else {
                    value = conf_p[conf_from, conf_to]
                }

                if (c != 'prederr' && !(conf_from %in% rownames(conf_p) && conf_to %in% colnames(conf_p))) {
                    next
                } else if (c == 'prederr' && !(conf_to %in% names(conf_p))) {
                    next
                }

                edge_info = edge_label(value, c, e, mod, p.min)
                eAttrs$label[e] = edge_info$label
                eAttrs$lwd[e] = edge_info$lwd
                eAttrs$fontcolor[e] = edge_info$fontcolor
            }
        }        
    }
    cat('RGraphviz object prepared.\n')

    ## Remove arrows from logic node (hidden and).
    
    for (e in edge_names) {
        edge = unlist(strsplit(e, '~'))
        from = edge[1]
        to = edge[2]

        if (is.logic.node.down(to)) {
            eAttrs$logic[e] = TRUE
            eAttrs$arrowsize[e] = 0

            if (substr(to, start = 1, stop = 2) == 'OR')
                eAttrs$color[e] = 'orange'
            if (substr(to, start = 1, stop = 3) == 'XOR')
                eAttrs$color[e] = 'red'
            if (substr(to, start = 1, stop = 3) == 'AND')
                eAttrs$color[e] = 'darkgreen'

            eAttrs$lty[e] = 'dashed'

            nAttrs$shape[to] = 'circle' 
        }

        if (is.logic.node.up(from)) {
            eAttrs$logic[e] = TRUE
            eAttrs$arrowsize[e] = 0

            eAttrs$lty[e] = 'dashed'

            if (substr(from, start = 1, stop = 4) == 'UPOR')
                eAttrs$color[e] = 'orange'
            if (substr(from, start = 1, stop = 5) == 'UPXOR')
                eAttrs$color[e] = 'red'
            if (substr(from, start = 1, stop = 5) == 'UPAND')
                eAttrs$color[e] = 'darkgreen'


        } else if (substr(from, start = 1, stop = 1) == '*') {
            eAttrs$logic[e] = TRUE
            eAttrs$arrowsize[e] = 0
            eAttrs$color[e] = 'black'
        }
    }

    if (pf) {
        cat('*** Add prima facie edges: ')
        ## For each edge...
        
        bic = adj.matrix$adj.matrix.bic

        for(e in edge_names) {
            edge = unlist(strsplit(e, '~'))
            from = edge[1]
            old_name = hypos_new_name[[from]]
            if (!is.null(old_name)) {
                from = old_name
            }
            to = edge[2]
            if (substr(to, start = 1, stop = 1) == '*') {
                to = substr(to, start = 3, stop = nchar(to))
            }

            ## ...checks if edge is present in BIC
            ## Check if edge in BIC (valid only if not logic edge) and
            ## 'to' is not a fake and
            
            if ((from %in% rownames(bic))
                && (to %in% colnames(bic))
                && !eAttrs$logic[e]
                && bic[from, to] == 0) {
                eAttrs$color[e] = 'red'
            } else {
                ## No PF
            }
        }
        cat('done')
    }

    if (sec) {
        pri.adj = primary$adj.matrix$adj.matrix.fit
        for (from in rownames(pri.adj)) {
            for (to in colnames(pri.adj)) {
                from.alt.name = from
                to.alt.name = to
                if (from %in% hypos_new_name) {
                    matching_nodes = names(which(hypos_new_name == from))
                    for (node in matching_nodes) {
                        if (is.logic.node.down(node)) {
                            from.alt.name = node
                        }
                    }
                }

                if (to %in% hypos_new_name) {
                    matching_nodes = names(which(hypos_new_name == to))
                    for(node in matching_nodes) {
                        if (is.logic.node.up(node)) {
                            to.alt.name = node
                        }
                    }
                }

                if (pri.adj[from, to] == 1) {
                    eAttrs$color[paste(from.alt.name, to.alt.name, sep='~')] = edge.color
                }

            }
        }
    }

    cat('Plotting graph and adding legends.\n')
    plot(graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, main = title, ... )

    ## Adds the legend to the plot.
    
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

        ## Add thickness legend.
        valid_names = node_names

        if (!disconnected) { 
            del = which(rowSums(hypo_mat) + colSums(hypo_mat) == 0 )
            w = !(rownames(hypo_mat) %in% names(del))
            valid_names = rownames(hypo_mat[w,])
        }

        if (expand) {
            valid_names =
                valid_names[unlist(lapply(valid_names,
                                          function(x) {
                                              !is.logic.node(x)
                                          }))]
        }

        valid_names = grep('^[*]_(.+)$', valid_names, value = TRUE, invert=TRUE)


        text = ""
        stat.pch = 0
        pt.bg = "white"
        col = "white"
        eloss = FALSE
        if ('eloss' %in% confidence) {
            eloss = TRUE
        }
        confidence = confidence[confidence != 'eloss']
        if (any(!is.na(confidence))) {
            text =
                c(expression(bold('Edge confidence')),
                  lapply(confidence,
                         function(x) {
                             if (x == "hg")
                                 return("Hypergeometric test")
                             if (x == "tp")
                                 return("Temporal Priority")
                             if (x == "pr")
                                 return("Probability Raising")
                             if (x == "pb")
                                 return("Parametric Bootstrap")
                             if (x == "sb")
                                 return("Statistical Bootstrap")
                             if (x == "npb")
                                 return("Non Parametric Bootstrap")
                             if (x == "prederr")
                                 return("Prediction Error")
                             if (x == "posterr")
                                 return("Posterior Classification Error")

                         }),
                  paste("p-value cutoff <", p.min))
            
            stat.pch = c(0, rep(18, length(confidence)), 0)
            pt.bg = c('white', rep('white', length(confidence)), 'white')
            col = c('white', rep('black', length(confidence)), 'white')
        }

        if ('Pattern' %in% as.types(x)) {
            y = delete.model(x)
            y = delete.type(y, 'Pattern')
        } else 
            y = x

        text =
            c(text, 
              ' ',
              expression(bold('Sample size')),
              paste0('n = ', nsamples(x), ', m = ', nevents(x)),
              paste0('|G| = ', ngenes(y), ', |P| = ', npatterns(x))     
              ) 

        stat.pch = c(stat.pch, rep(0, 2), rep(20, 2), rep(0, 2))
        pt.bg = c(pt.bg, rep('white', 2), rep('black', 2), rep('white', 2))
        col = c(col, rep('white', 2), rep('white', 2), rep('white', 2))

        mods = NULL

        if (eloss) {
            for (model in models) {
                mods_label = gsub('_', ' ', model)
                if (!is.null(x$kfold) && !is.null(get(model, x$kfold)$eloss)) {
                    mods_label = paste(mods_label, '- eloss:', round(mean(get(model, x$kfold)$eloss), 5))
                }
                mods = c(mods, mods_label)
            }
        } else {
            for (model in models) {
                mods_label = gsub('_', ' ', model)
                mods = c(mods, mods_label)
            }
        }
        text = 
            c(text, '\n',
              expression(bold('Algorithm:')),
              paste0(mods)
              )

        stat.pch = c(stat.pch, 20)
        pt.bg = c(pt.bg, 'black')
        col = c(col,'black') 

        if (length(models) > 1) {
            stat.pch = c(stat.pch, 20)
            pt.bg = c(pt.bg, 'darkgrey')
            col = c(col, 'darkgrey') 
        }


        legend(legend.pos.l,
               legend = text,
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

    if (!is.na(file)) {
        cat('Saving visualized device to file:', file)
        dev.copy2pdf(file = file)
    }
    cat('\n')

    if (export.igraph) {
        output = list()
        output$graph = igraph.from.graphNEL(graph)
        output$nodes = nAttrs
        output$edges = eAttrs
        output$description = title
        output$models = models
        return(output)
    }
}

#### end of file -- tronco.R
