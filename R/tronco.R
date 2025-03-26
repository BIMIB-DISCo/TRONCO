#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2017, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Reconstruct a progression model using CAPRESE algorithm. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.caprese(test_dataset_no_hypos)
#'
#' @title tronco caprese
#' @param data A TRONCO compliant dataset.
#' @param lambda Coefficient to combine the raw estimate with a correction factor into a shrinkage estimator.
#' @param silent A parameter to disable/enable verbose messages.
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.caprese
#' @importFrom stats phyper
#' @importFrom bnlearn empty.graph set.arc
#' @importFrom igraph graph.adjacency get.shortest.paths
#'
tronco.caprese <- function(data,
                           lambda = 0.5,
                           silent = FALSE,
                           epos = 0.0,
                           eneg = 0.0) {
  ## Check for the inputs to be correct
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
  }
  if (lambda < 0 || lambda > 1) {
    stop("The value of the shrinkage parameter lambda has to be in [0:1]!",
         call. = FALSE)
    
  }
  
  ## check for the input to be compliant
  is.compliant(data)
  
  ## check if there are hypotheses
  if (npatterns(data) > 0) {
    warning("Patters found in input for tronco.caprese\n")
  }

  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  ## Reconstruct the reconstruction with CAPRESE
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: CAPRESE with shrinkage coefficient: ',
        lambda,
        '.\n'
      )
    )
  }
  reconstruction = caprese.fit(
    dataset = data$genotypes,
    lambda = lambda,
    silent = silent,
    epos = epos,
    eneg = eneg,
    hypotheses = data$hypotheses
  )
  
  ## Structure to save the results
  results = data
  
  results$confidence = reconstruction$confidence
  
  results$model = reconstruction$model
  
  results$parameters = reconstruction$parameters
  
  results$execution.time = reconstruction$execution.time
  
  
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
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
  
}


#' Reconstruct a progression model using CAPRI algorithm. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset)
#' recon = tronco.capri(test_dataset, nboot = 1)
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
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @param restart An integer, the number of random restarts.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.capri
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths is.dag
#' @importFrom stats phyper AIC BIC wilcox.test
#'
tronco.capri <- function(data,
                         command = "hc",
                         regularization = c("bic", "aic"),
                         do.boot = TRUE,
                         nboot = 100,
                         pvalue = 0.05,
                         min.boot = 3,
                         min.stat = TRUE,
                         boot.seed = NULL,
                         silent = FALSE,
                         epos = 0.0,
                         eneg = 0.0,
                         restart = 100) {
  ## Check for the inputs to be correct
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
    
  }
  
  ## Enforce data to be numeric
  data = enforce.numeric(data)
  
  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  if (!command %in% c("hc", "tabu")) {
    stop("The inference can be performed either by hill climbing or tabu search!",
         call. = FALSE)
  }
  
  if (pvalue < 0 || pvalue > 1) {
    stop("The value of the pvalue has to be in [0:1]!", call. = FALSE)
    
  }
  
  if (!all(regularization %in% c('loglik', 'bic', 'aic'))) {
    stop("Possible regularization are loglik, bic or aic", call. = FALSE)
    
  }
  
  if (epos < 0 || epos >= 0.5 || eneg < 0 || eneg >= 0.5) {
    stop("The values of the error rates have to be in [0:0.5)!",
         call. = FALSE)
  }
  
  ## Check for the input to be compliant
  is.compliant(data)
  
  ## Reconstruct the reconstruction with CAPRI
  if (is.null(boot.seed)) {
    my.seed = "NULL"
  }
  else {
    my.seed = boot.seed
    
  }
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: CAPRI with \"',
        paste0(regularization, collapse = ", "),
        '\" regularization and \"',
        command,
        '\" likelihood-fit strategy.\n',
        '\tRandom seed: ',
        my.seed,
        '.\n',
        '\tBootstrap iterations (Wilcoxon): ',
        ifelse(do.boot, nboot, 'disabled'),
        '.\n',
        ifelse(
          do.boot,
          paste0(
            '\t\texhaustive bootstrap: ',
            min.stat,
            '.\n\t\tp-value: ',
            pvalue,
            '.\n\t\tminimum bootstrapped scores: ',
            min.boot,
            '.\n'
          ),
          ''
        )
      )
    )
  }
  
  reconstruction =
    capri.fit(
      data$genotypes,
      data$hypotheses,
      command = command,
      regularization = regularization,
      do.boot = do.boot,
      nboot = nboot,
      pvalue = pvalue,
      min.boot = min.boot,
      min.stat = min.stat,
      boot.seed = boot.seed,
      silent = silent,
      epos = epos,
      eneg = eneg,
      restart = restart
    )
  
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
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
  
}


#' Reconstruct a progression model using Edmonds algorithm combined
#' with probabilistic causation. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.edmonds(test_dataset_no_hypos, nboot = 1)
#'
#' @title Tronco Edmonds
#' @param data A TRONCO compliant dataset.
#' @param regularization Select the regularization for the
#' likelihood estimation, e.g., BIC, AIC.
#' @param score Select the score for the estimation of
#' the best tree, e.g., pointwise mutual information (pmi), conditional entropy (entropy).
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
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.edmonds
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths is.dag
#### @importFrom infotheo mutinformation
#' @importFrom stats phyper AIC BIC
#'
tronco.edmonds <- function(data,
                           regularization = "no_reg",
                           score = "pmi",
                           do.boot = TRUE,
                           nboot = 100,
                           pvalue = 0.05,
                           min.boot = 3,
                           min.stat = TRUE,
                           boot.seed = NULL,
                           silent = FALSE,
                           epos = 0.0,
                           eneg = 0.0) {
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
    
  }
  
  ## Enforce data to be numeric
  data = enforce.numeric(data)
  
  ## Check for the inputs to be correct.
  
  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  if (pvalue < 0 || pvalue > 1) {
    stop("The value of the pvalue has to be in [0:1]!", call. = FALSE)
    
  }
  
  if (!all(regularization %in% c('no_reg', 'loglik', 'bic', 'aic'))) {
    stop("Possible regularization are no-reg, loglik, bic or aic",
         call. = FALSE)
    
  }
  
  if (!all(score %in% c('pmi', 'mi', 'entropy', 'cpmi'))) {
    stop("Possible scores are pmi, mi, entropy or cpmi", call. = FALSE)
    
  }
  
  if (epos < 0 || epos >= 0.5 || eneg < 0 || eneg >= 0.5) {
    stop("The values of the error rates have to be in [0:0.5)!",
         call. = FALSE)
  }
  
  ## Check for the input to be compliant.
  
  is.compliant(data)
  
  ## check if there are hypotheses
  
  if (npatterns(data) > 0) {
    warning("Patters found in input for tronco.edmonds\n")
  }
  
  ## Reconstruct the reconstruction with Edmonds.
  
  if (is.null(boot.seed)) {
    my.seed = "NULL"
  }
  else {
    my.seed = boot.seed
    
  }
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: Edmonds with \"',
        paste0(regularization, collapse = ", "),
        '\" regularization',
        '\tRandom seed: ',
        my.seed,
        '.\n',
        '\tBootstrap iterations (Wilcoxon): ',
        ifelse(do.boot, nboot, 'disabled'),
        '.\n',
        ifelse(
          do.boot,
          paste0(
            '\t\texhaustive bootstrap: ',
            min.stat,
            '.\n\t\tp-value: ',
            pvalue,
            '.\n\t\tminimum bootstrapped scores: ',
            min.boot,
            '.\n'
          ),
          ''
        )
      )
    )
  }
  
  reconstruction =
    edmonds.fit(
      data$genotypes,
      regularization = regularization,
      score = score,
      do.boot = do.boot,
      nboot = nboot,
      pvalue = pvalue,
      min.boot = min.boot,
      min.stat = min.stat,
      boot.seed = boot.seed,
      silent = silent,
      epos = epos,
      eneg = eneg,
      hypotheses = data$hypotheses
    )
  
  ## Structure to save the results.
  results = data
  results$adj.matrix.prima.facie = reconstruction$adj.matrix.prima.facie
  results$adj.matrix.prima.facie.cyclic = reconstruction$adj.matrix.prima.facie.cyclic
  results$confidence = reconstruction$confidence
  results$model = reconstruction$model
  results$parameters = reconstruction$parameters
  results$execution.time = reconstruction$execution.time
  
  ## Add BIC/AIC/LogLik informations
  
  if (!silent) {
    cat('*** Evaluating BIC / AIC / LogLik informations.\n')
  }
  
  search_scores = score
  
  if ("no_reg" %in% regularization) {
    for (my_s in search_scores) {
      bayes.net = as.bnlearn.network(results, model = paste('edmonds_no_reg', my_s, sep =
                                                              "_"))
      score = logLik(bayes.net$net, data = bayes.net$data)
      logLik = score
      results$model[[paste('edmonds_no_reg', my_s, sep = "_")]]$score = score
      results$model[[paste('edmonds_no_reg', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("loglik" %in% regularization) {
    for (my_s in search_scores) {
      bayes.net = as.bnlearn.network(results, model  = paste('edmonds_loglik', my_s, sep =
                                                               "_"))
      score = logLik(bayes.net$net, data = bayes.net$data)
      logLik = score
      results$model[[paste('edmonds_loglik', my_s, sep = "_")]]$score = score
      results$model[[paste('edmonds_loglik', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("bic" %in% regularization) {
    for (my_s in search_scores) {
      bayes.net = as.bnlearn.network(results, model  = paste('edmonds_bic', my_s, sep =
                                                               "_"))
      score = BIC(bayes.net$net, data = bayes.net$data)
      logLik = logLik(bayes.net$net, data = bayes.net$data)
      results$model[[paste('edmonds_bic', my_s, sep = "_")]]$score = score
      results$model[[paste('edmonds_bic', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("aic" %in% regularization) {
    for (my_s in search_scores) {
      bayes.net = as.bnlearn.network(results, model = paste('edmonds_aic', my_s, sep =
                                                              "_"))
      score = AIC(bayes.net$net, data = bayes.net$data)
      logLik = logLik(bayes.net$net, data = bayes.net$data)
      results$model[[paste('edmonds_aic', my_s, sep = "_")]]$score = score
      results$model[[paste('edmonds_aic', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  ## the reconstruction has been completed
  if (!silent)
    cat(paste(
      "The reconstruction has been successfully completed in",
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
  
}


#' Reconstruct a progression model using Gabow algorithm combined
#' with probabilistic causation. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.gabow(test_dataset_no_hypos, nboot = 1)
#'
#' @title Tronco Gabow
#' @param data A TRONCO compliant dataset.
#' @param regularization Select the regularization for the
#' likelihood estimation, e.g., BIC, AIC.
#' @param score Select the score for the estimation of
#' the best tree, e.g., pointwise mutual information (pmi), conditional entropy (entropy).
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
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @param do.raising Whether to use or not the raising condition as a prior.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.gabow
#' @importFrom bnlearn hc tabu empty.graph set.arc score amat<- amat
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths graph_from_adjacency_matrix clusters unfold.tree
#' @importFrom igraph is.dag
#' @importFrom gtools permutations
#' @importFrom stats phyper AIC BIC logLik runif
#'
tronco.gabow <- function(data,
                         regularization = "no_reg",
                         score = "pmi",
                         do.boot = TRUE,
                         nboot = 100,
                         pvalue = 0.05,
                         min.boot = 3,
                         min.stat = TRUE,
                         boot.seed = NULL,
                         silent = FALSE,
                         epos = 0.0,
                         eneg = 0.0,
                         do.raising = TRUE) {
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
    
  }
  
  ## Enforce data to be numeric
  data = enforce.numeric(data)
  
  ## Check for the inputs to be correct.
  
  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  if (pvalue < 0 || pvalue > 1) {
    stop("The value of the pvalue has to be in [0:1]!", call. = FALSE)
    
  }
  
  if (!all(regularization %in% c('no_reg', 'loglik', 'bic', 'aic'))) {
    stop("Possible regularization are no-reg, loglik, bic or aic",
         call. = FALSE)
    
  }
  
  if (!all(score %in% c('pmi', 'mi', 'entropy', 'cpmi'))) {
    stop("Possible scores are pmi, mi, entropy or cpmi", call. = FALSE)
  }
  
  if (epos < 0 || epos >= 0.5 || eneg < 0 || eneg >= 0.5) {
    stop("The values of the error rates have to be in [0:0.5)!",
         call. = FALSE)
  }
  
  ## Check for the input to be compliant.
  
  is.compliant(data)
  
  ## check if there are hypotheses
  
  if (npatterns(data) > 0) {
    warning("Patters found in input for tronco.gabow\n")
  }
  
  ## Reconstruct the reconstruction with MLE.
  
  if (is.null(boot.seed)) {
    my.seed = "NULL"
  }
  else {
    my.seed = boot.seed
    
  }
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: Gabow with \"',
        paste0(regularization, collapse = ", "),
        '\" regularization',
        '\tRandom seed: ',
        my.seed,
        '.\n',
        '\tBootstrap iterations (Wilcoxon): ',
        ifelse(do.boot, nboot, 'disabled'),
        '.\n',
        ifelse(
          do.boot,
          paste0(
            '\t\texhaustive bootstrap: ',
            min.stat,
            '.\n\t\tp-value: ',
            pvalue,
            '.\n\t\tminimum bootstrapped scores: ',
            min.boot,
            '.\n'
          ),
          ''
        )
      )
    )
  }
  
  reconstruction =
    gabow.fit(
      data$genotypes,
      regularization = regularization,
      score = score,
      do.boot = do.boot,
      nboot = nboot,
      pvalue = pvalue,
      min.boot = min.boot,
      min.stat = min.stat,
      boot.seed = boot.seed,
      silent = silent,
      epos = epos,
      eneg = eneg,
      do.raising = do.raising,
      hypotheses = data$hypotheses
    )
  
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
  
  search_scores = score
  
  is.acyclic = TRUE
  
  models.adj.matrix = as.adj.matrix(results)
  
  if ("no_reg" %in% regularization) {
    for (my_s in search_scores) {
      mod.name = paste('gabow_no_reg', my_s, sep = "_")
      this.matrix = models.adj.matrix[[mod.name]]
      if (is.dag(graph.adjacency(this.matrix))) {
        bayes.net = as.bnlearn.network(results, model = mod.name)
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
      } else {
        score = -1
        logLik = -1
        is.acyclic = FALSE
      }
      results$model[[paste('gabow_no_reg', my_s, sep = "_")]]$score = score
      results$model[[paste('gabow_no_reg', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("loglik" %in% regularization) {
    for (my_s in search_scores) {
      mod.name = paste('gabow_loglik', my_s, sep = "_")
      this.matrix = models.adj.matrix[[mod.name]]
      if (is.dag(graph.adjacency(this.matrix))) {
        bayes.net = as.bnlearn.network(results, model = mod.name)
        score = logLik(bayes.net$net, data = bayes.net$data)
        logLik = score
      } else {
        score = -1
        logLik = -1
        is.acyclic = FALSE
      }
      results$model[[paste('gabow_loglik', my_s, sep = "_")]]$score = score
      results$model[[paste('gabow_loglik', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("bic" %in% regularization) {
    for (my_s in search_scores) {
      mod.name = paste('gabow_bic', my_s, sep = "_")
      this.matrix = models.adj.matrix[[mod.name]]
      if (is.dag(graph.adjacency(this.matrix))) {
        bayes.net = as.bnlearn.network(results, model = mod.name)
        score = BIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
      } else {
        score = -1
        logLik = -1
        is.acyclic = FALSE
      }
      results$model[[paste('gabow_bic', my_s, sep = "_")]]$score = score
      results$model[[paste('gabow_bic', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  if ("aic" %in% regularization) {
    for (my_s in search_scores) {
      mod.name = paste('gabow_aic', my_s, sep = "_")
      this.matrix = models.adj.matrix[[mod.name]]
      
      if (is.dag(graph.adjacency(this.matrix))) {
        bayes.net = as.bnlearn.network(results, model = mod.name)
        score = AIC(bayes.net$net, data = bayes.net$data)
        logLik = logLik(bayes.net$net, data = bayes.net$data)
      } else {
        score = -1
        logLik = -1
        is.acyclic = FALSE
      }
      results$model[[paste('gabow_aic', my_s, sep = "_")]]$score = score
      results$model[[paste('gabow_aic', my_s, sep = "_")]]$logLik = logLik
    }
  }
  
  ## the reconstruction has been completed.
  
  if (!is.acyclic) {
    save(results, file = paste0('result_',
                                as.character(as.integer(runif(
                                  1
                                ) * 10000)),
                                '.RData'))
  }
  
  if (!silent)
    cat(paste(
      "The reconstruction has been successfully completed in",
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
}


#' Reconstruct a progression model using Chow Liu
#' algorithm combined with probabilistic causation. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.chowliu(test_dataset_no_hypos, nboot = 1)
#'
#' @title Tronco Chow Liu
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
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @return A TRONCO compliant object with reconstructed
#' model
#' @export tronco.chowliu
#' @importFrom bnlearn hc tabu empty.graph set.arc chow.liu
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths is.dag
#' @importFrom stats phyper AIC BIC
#'
tronco.chowliu <- function(data,
                           regularization = c("bic", "aic"),
                           do.boot = TRUE,
                           nboot = 100,
                           pvalue = 0.05,
                           min.boot = 3,
                           min.stat = TRUE,
                           boot.seed = NULL,
                           silent = FALSE,
                           epos = 0.0,
                           eneg = 0.0) {
  ## Check for the inputs to be correct.
  
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
    
  }
  
  ## Enforce data to be numeric
  data = enforce.numeric(data)
  
  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  if (pvalue < 0 || pvalue > 1) {
    stop("The value of the pvalue has to be in [0:1]!", call. = FALSE)
    
  }
  
  if (!all(regularization %in% c('loglik', 'bic', 'aic'))) {
    stop("Possible regularization are loglik, bic or aic", call. = FALSE)
    
  }
  
  if (epos < 0 || epos >= 0.5 || eneg < 0 || eneg >= 0.5) {
    stop("The values of the error rates have to be in [0:0.5)!",
         call. = FALSE)
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
    my.seed = boot.seed
    
  }
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: Chow Liu with \"',
        paste0(regularization, collapse = ", "),
        '\" regularization',
        '\tRandom seed: ',
        my.seed,
        '.\n',
        '\tBootstrap iterations (Wilcoxon): ',
        ifelse(do.boot, nboot, 'disabled'),
        '.\n',
        ifelse(
          do.boot,
          paste0(
            '\t\texhaustive bootstrap: ',
            min.stat,
            '.\n\t\tp-value: ',
            pvalue,
            '.\n\t\tminimum bootstrapped scores: ',
            min.boot,
            '.\n'
          ),
          ''
        )
      )
    )
  }
  
  reconstruction =
    chow.liu.fit(
      data$genotypes,
      regularization = regularization,
      do.boot = do.boot,
      nboot = nboot,
      pvalue = pvalue,
      min.boot = min.boot,
      min.stat = min.stat,
      boot.seed = boot.seed,
      silent = silent,
      epos = epos,
      eneg = eneg,
      hypotheses = data$hypotheses
    )
  
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
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
  
}


#' Reconstruct a progression model using Prim algorithm combined with probabilistic causation. For details and examples
#' regarding the inference process and on the algorithm implemented in the package,
#' we refer to the Vignette Section 6.
#'
#' @examples
#' data(test_dataset_no_hypos)
#' recon = tronco.prim(test_dataset_no_hypos, nboot = 1)
#'
#' @title Tronco Prim
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
#' @param epos Error rate of false positive errors.
#' @param eneg Error rate of false negative errors.
#' @return A TRONCO compliant object with reconstructed model
#' @export tronco.prim
#' @importFrom bnlearn hc tabu empty.graph set.arc
#' @importFrom igraph get.edgelist E E<-
#' @importFrom igraph graph.adjacency get.adjacency graph.union edge
#' @importFrom igraph get.shortest.paths minimum.spanning.tree is.dag
#### @importFrom infotheo mutinformation
#' @importFrom stats phyper AIC BIC
#'
tronco.prim <- function(data,
                        regularization = "no_reg",
                        do.boot = TRUE,
                        nboot = 100,
                        pvalue = 0.05,
                        min.boot = 3,
                        min.stat = TRUE,
                        boot.seed = NULL,
                        silent = FALSE,
                        epos = 0.0,
                        eneg = 0.0) {
  ## Check for the inputs to be correct.
  
  if (is.null(data) || is.null(data$genotypes)) {
    stop("The dataset given as input is not valid.")
    
  }
  
  ## Enforce data to be numeric
  data = enforce.numeric(data)
  
  if (is.null(data$hypotheses)) {
    data$hypotheses = NA
    
  }
  
  if (pvalue < 0 || pvalue > 1) {
    stop("The value of the pvalue has to be in [0:1]!", call. = FALSE)
    
  }
  
  if (!all(regularization %in% c('no_reg', 'loglik', 'bic', 'aic'))) {
    stop("Possible regularization are no-reg, loglik, bic or aic",
         call. = FALSE)
    
  }
  
  if (epos < 0 || epos >= 0.5 || eneg < 0 || eneg >= 0.5) {
    stop("The values of the error rates have to be in [0:0.5)!",
         call. = FALSE)
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
    my.seed = boot.seed
    
  }
  if (silent == FALSE) {
    cat('*** Checking input events.\n')
    invalid = consolidate.data(data, TRUE)
    if (length(unlist(invalid)) > 0)
      warning("Input events should be consolidated - see consolidate.data.")
    
    
    cat(
      paste0(
        '*** Inferring a progression model with the following settings.\n',
        '\tDataset size: n = ',
        nsamples(data),
        ', m = ',
        nevents(data),
        '.\n',
        '\tAlgorithm: Prim with \"',
        paste0(regularization, collapse = ", "),
        '\" regularization',
        '\tRandom seed: ',
        my.seed,
        '.\n',
        '\tBootstrap iterations (Wilcoxon): ',
        ifelse(do.boot, nboot, 'disabled'),
        '.\n',
        ifelse(
          do.boot,
          paste0(
            '\t\texhaustive bootstrap: ',
            min.stat,
            '.\n\t\tp-value: ',
            pvalue,
            '.\n\t\tminimum bootstrapped scores: ',
            min.boot,
            '.\n'
          ),
          ''
        )
      )
    )
  }
  
  reconstruction =
    prim.fit(
      data$genotypes,
      regularization = regularization,
      do.boot = do.boot,
      nboot = nboot,
      pvalue = pvalue,
      min.boot = min.boot,
      min.stat = min.stat,
      boot.seed = boot.seed,
      silent = silent,
      epos = epos,
      eneg = eneg,
      hypotheses = data$hypotheses
    )
  
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
    bayes.net = as.bnlearn.network(results, model = 'prim_no_reg')
    score = logLik(bayes.net$net, data = bayes.net$data)
    logLik = score
    results$model$prim_no_reg$score = score
    results$model$prim_no_reg$logLik = logLik
  }
  
  if ("loglik" %in% regularization) {
    bayes.net = as.bnlearn.network(results, model = 'prim_loglik')
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
  
  ## the reconstruction has been completed
  if (!silent)
    cat(paste(
      "The reconstruction has been successfully completed in",
      format(.POSIXct(
        round(reconstruction$execution.time[3],
              digits = 0),
        tz = "GMT"
      ),
      "%Hh:%Mm:%Ss"),
      "\n"
    ))
  
  
  return(results)
  
}


#' Bootstrap a reconstructed progression model. For details and examples
#' regarding the statistical assesment of an inferred model,
#' we refer to the Vignette Section 7.
#'
#' @examples
#' data(test_model)
#' boot = tronco.bootstrap(test_model, nboot = 1, cores.ratio = 0)
#'
#' @title tronco bootstrap
#' @param reconstruction The output of tronco.capri or
#' tronco.caprese
#' @param type Parameter to define the type of sampling
#' to be performed, e.g., non-parametric for uniform sampling.
#' @param nboot Number of bootstrap sampling to be performed
#' when estimating the model confidence.
#' @param cores.ratio Percentage of cores to use
#' coresRate * (numCores - 1)
#' @param silent A parameter to disable/enable verbose messages.
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
                             cores.ratio = 1,
                             silent = FALSE) {
  ## Check for the input to be compliant.
  is.compliant(reconstruction)
  
  ## Check for the inputs to be given.
  is.model(reconstruction)
  
  
  if (type == "statistical"
      && !((
        reconstruction$parameters$algorithm == "CAPRI"
        || reconstruction$parameters$algorithm == "PRIM"
        || reconstruction$parameters$algorithm == "CHOW_LIU"
        || reconstruction$parameters$algorithm == "EDMONDS"
      )
      && reconstruction$parameters$do.boot == TRUE
      )) {
    stop(
      paste(
        "To perform statistical bootstrap, the algorithm used for",
        "the reconstruction must be CAPRI, PRIM, CHOW_LIU or EDMONDS",
        "with bootstrap."
      ),
      call. = FALSE
    )
  }
  
  ## Set all the needed parameters to perform the bootstrap estimation
  if (!type %in% c("non-parametric", "statistical")) {
    stop(
      paste(
        "The types of bootstrap that can be performed are:",
        "non-parametric or statistical."
      ),
      call. = FALSE
    )
  }
  
  ## Perform the selected bootstrap procedure
  if (!silent) {
    cat("*** Executing now the bootstrap procedure, this may take a long time...\n")
  }
  
  parameters = as.parameters(reconstruction)
  
  if (parameters$algorithm == "CAPRESE") {
    lambda = parameters$lambda
    curr.boot = bootstrap(reconstruction,
                          type,
                          nboot,
                          cores.ratio,
                          silent = silent)
    if (!silent) {
      cat(
        "Performed",
        type,
        "bootstrap with",
        nboot,
        "resampling and",
        lambda,
        "as shrinkage parameter.\n"
      )
    }
    
  } else {
    curr.boot = bootstrap(reconstruction,
                          type,
                          nboot,
                          cores.ratio,
                          silent = silent)
    if (!silent) {
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
  }
  reconstruction$bootstrap = curr.boot
  return(reconstruction)
}


#### end of file -- tronco.R
