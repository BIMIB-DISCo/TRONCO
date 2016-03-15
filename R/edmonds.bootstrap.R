#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# perform non-parametric or parametric bootstrap to evalutate the confidence of the reconstruction
# @title bootstrap.edmonds
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot.algorithm integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param do.estimation should I perform the estimation of the error rates and probabilities?
# @param silent should I be verbose?
# @param reconstruction Result of a previous reconstruction
# @param command should I perform non-parametric or parametric bootstrap?
# @param nboot number of bootstrap resampling to be performed
# @param bootstrap.statistics Result of a previous bootstrap analysis
# @param verbose Should I print messages?
# @return bootstrap.statistics: statistics of the bootstrap
#
bootstrap.edmonds <- function(dataset, 
                              regularization, 
                              do.boot,
                              nboot.algorithm, 
                              pvalue,
                              min.boot,
                              min.stat,
                              boot.seed,
                              do.estimation,
                              silent,
                              reconstruction, 
                              command = "non-parametric",
                              nboot = 100,
                              bootstrap.statistics = list(),
                              verbose = FALSE,
                              cores.ratio = 1) {
    
    ## Start the clock to measure the execution time
    
    ## library(doParallel)

    ptm <- proc.time();
    
    ## Structure to save the results of the bootstrap.
    
    curr.bootstrap.results = array(list(-1), c(nboot,nevents(reconstruction)))
    colnames(curr.bootstrap.results) = rownames(as.events(reconstruction))
    bootstrap.results = list()
    bootstrap.results[names(as.models(reconstruction))] = list(curr.bootstrap.results)
    
    curr.bootstrap.adj.matrix = array(list(0), c(nevents(reconstruction)+1,nevents(reconstruction)+1))
    colnames(curr.bootstrap.adj.matrix) = c("None",rownames(as.events(reconstruction)))
    rownames(curr.bootstrap.adj.matrix) = c("None",rownames(as.events(reconstruction)))
    bootstrap.adj.matrix = list()
    bootstrap.adj.matrix[names(as.models(reconstruction))] = list(curr.bootstrap.adj.matrix)
    
    bootstrap.adj.matrix.frequency = list()
    bootstrap.adj.matrix.frequency[names(as.models(reconstruction))] = list(curr.bootstrap.adj.matrix)
    
    curr.edge.confidence = array(list(0), c(nevents(reconstruction),nevents(reconstruction)))
    colnames(curr.edge.confidence) = rownames(as.events(reconstruction))
    rownames(curr.edge.confidence) = rownames(as.events(reconstruction))
    bootstrap.edge.confidence = list()
    bootstrap.edge.confidence[names(as.models(reconstruction))] = list(curr.edge.confidence)
    
    overall.confidence = list()
    overall.confidence[names(as.models(reconstruction))] = list(0)
    overall.frequency = list()
    overall.frequency[names(as.models(reconstruction))] = list(0)

    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }



    expected.execution.time =
        round(((reconstruction$execution.time[3] * nboot) / (cores)), digits = 0)
    cat("Expected completion in approx.",
        format(.POSIXct(expected.execution.time, tz="GMT"),
               "%Hh:%Mm:%Ss"),
        "\n")

    if (!verbose) {
        cl = makeCluster(cores)    
    } else {
        cl = makeCluster(cores, outfile="")
    }

    registerDoParallel(cl)
    cat('*** Using', cores, 'cores via "parallel" \n')
    
    ## Perform nboot bootstrap resampling
    ## for (num in 1:nboot

    ## If parametric bootstrap selected prepare input.

    if (command == 'parametric') {
        
        ## Structure to save the samples probabilities.
        
        samples.probabilities = list();
        
        ## Define the possible samples given the current number of events.
        
        possible.strings = 2 ^ ncol(dataset)
        
        err = ""
        message = "Too many events in the dataset! Parametric bootstrastap can not be performed."
        err =
            tryCatch({
                curr.dataset = suppressWarnings(array(0, c(possible.strings, ncol(dataset))))
            }, error = function(e) {
                err <- message
            })
        
        
        if (toString(err) == message) {
            stop(err, call. = FALSE)
        }
        
        for (i in 1:possible.strings) {
            curr.dataset[i, ] = decimal.to.binary.tree(i - 1, ncol(dataset))
        }

        colnames(curr.dataset) = colnames(dataset)
        
        for (m in names(as.models(reconstruction))) {
            
            ## Estimate the samples probabilities for each model.
            
            samples.probabilities[m] =
                list(estimate.tree.samples(curr.dataset,
                                          as.adj.matrix(reconstruction, models = m)[[m]],
                                          as.marginal.probs(reconstruction, models = m, type = "fit")[[m]],
                                          as.conditional.probs(reconstruction, models = m, type = "fit")[[m]],
                                          as.error.rates(reconstruction, models = m)[[m]]))
        }
    }
    
    r =
        foreach(num = 1:nboot ) %dopar% {    
            
            ## Performed the bootstrapping procedure.
            
            if (command == "non-parametric") {
                
                ## Perform the sampling for the current step of bootstrap.
                
                samples = sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE)
                bootstrapped.dataset = dataset[samples,]
                
                curr.reconstruction = list()
                curr.reconstruction$genotypes = bootstrapped.dataset
                curr.reconstruction$annotations = reconstruction$annotations
                curr.reconstruction$types = reconstruction$types
                
                ## Perform the reconstruction on the bootstrapped dataset.
                
                bootstrapped.topology =
                    tronco.edmonds(curr.reconstruction,
                                 regularization,
                                 do.boot,
                                 nboot.algorithm,
                                 pvalue,
                                 min.boot,
                                 min.stat,
                                 boot.seed,
                                 silent)

                curr.reconstruction = bootstrapped.topology;

            } else if (command=="parametric") {
                
                ## Perform the reconstruction for each model.
                
                new.reconstruction = reconstruction;
                new.reconstruction$model = list();
                for (m in names(as.models(reconstruction))) {
                    ## Perform the sampling for the current step of
                    ## bootstrap and regularizator.
                    
                    samples =
                        sample(1:nrow(curr.dataset), 
                               size = nrow(dataset), 
                               replace = TRUE, 
                               prob = samples.probabilities[[m]])

                    bootstrapped.dataset = curr.dataset[samples,]
                    
                    curr.reconstruction = list()
                    curr.reconstruction$genotypes = bootstrapped.dataset
                    curr.reconstruction$annotations = reconstruction$annotations
                    curr.reconstruction$types = reconstruction$types
                    
                    ## Perform the reconstruction on the bootstrapped
                    ## dataset.
                    
                    bootstrapped.topology =
                        tronco.edmonds(curr.reconstruction,
                                     m,
                                     do.boot,
                                     nboot.algorithm,
                                     pvalue,
                                     min.boot,
                                     min.stat,
                                     boot.seed,
                                     silent)
                    
                    ## Save the results for this model.
                    
                    new.reconstruction$model[m] =
                        as.models(bootstrapped.topology,models = m)
                }
                curr.reconstruction = new.reconstruction;
            } else if (command == "statistical") {
                
                curr.reconstruction = list()
                curr.reconstruction$genotypes = reconstruction$genotypes;
                curr.reconstruction$annotations = reconstruction$annotations;
                curr.reconstruction$types = reconstruction$types;
                
                ## Perform the reconstruction on the bootstrapped
                ## dataset.
                
                bootstrapped.topology =
                    tronco.edmonds(curr.reconstruction,
                                 regularization,
                                 do.boot,
                                 nboot.algorithm,
                                 pvalue,
                                 min.boot,
                                 min.stat,
                                 boot.seed,
                                 silent)
                
                curr.reconstruction = bootstrapped.topology;
            }
            
            ## Set the reconstructed selective advantage edges.
            
            bootstrap.results = list()
            for (m in names(as.models(curr.reconstruction))) {
                
                ## Get the parents pos.
                parents.pos =
                    array(list(), c(nevents(curr.reconstruction), 1))
                
                
                curr.adj.matrix =
                    as.adj.matrix(curr.reconstruction, models = m)[[m]]
                for (i in 1:nevents(curr.reconstruction)) {
                    for (j in 1:nevents(curr.reconstruction)) {
                        if (i != j && curr.adj.matrix[i, j] == 1) {
                            parents.pos[j, 1] =
                                list(c(unlist(parents.pos[j, 1]), i))
                        }
                    }
                }
                
                parents.pos[unlist(lapply(parents.pos,is.null))] = list(-1)
                
                ## Save the results.
                
                bootstrap.results[[m]] = t(parents.pos);
            }
            bootstrap.results
        }

    stopCluster(cl)
    cat("\n*** Reducing results\n")

    for (m in names(bootstrap.results)) {
        y = Reduce(rbind, lapply(r, function(z, type) { get(type, z) }, type = m))
        bootstrap.results[[m]] = y
    }

    ## Set the statistics of the bootstrap.
    
    for (m in names(as.models(reconstruction))) {
        
        curr.bootstrap.adj.matrix = bootstrap.adj.matrix[[m]]
        
        for (i in 2:ncol(curr.bootstrap.adj.matrix)) {
            
            curr.result = bootstrap.results[[m]][ , i - 1]
            
            for (j in 1:length(curr.result)) {
                
                curr.val = curr.result[[j]]
                
                for (k in 1:length(curr.val)) {
                    if (length(curr.val[k]) == 1 && curr.val[k] == -1) {
                        curr.bootstrap.adj.matrix[[1, i]] =
                            curr.bootstrap.adj.matrix[[1, i]] + 1
                    } else {
                        curr.bootstrap.adj.matrix[[curr.val[k] + 1, i]] =
                            curr.bootstrap.adj.matrix[[curr.val[k] + 1, i]] + 1
                    }
                }
            }
            
        }
        
        bootstrap.adj.matrix[[m]] = curr.bootstrap.adj.matrix;
        rownames(bootstrap.results[[m]]) =
            paste("Iteration ", 1:nrow(bootstrap.results[[m]]), sep= "")
        
    }
    
    ## Evalutate the overall confidence.
    
    for (m in names(as.models(reconstruction))) {
        
        curr.bootstrap.results = bootstrap.results[[m]]
        
        for (i in 1:nrow(curr.bootstrap.results)) {
            
            curr.adj.matrix = array(0, c(nevents(reconstruction), nevents(reconstruction)))
            
            for (j in 1:ncol(curr.bootstrap.results)) {
                
                curr.result = curr.bootstrap.results[i, j]
                
                for (k in 1:length(curr.result)) {
                    
                    curr.val = curr.result[[k]]
                    
                    for (l in 1:length(curr.val)) {
                        if (length(curr.val[l]) > 1 || curr.val[l] != -1) {
                            curr.adj.matrix[curr.val[l], j] = 1
                        }
                    }
                }
            }
            
            ## If I have a perfect match between the reconstructed
            ## topologies, increase the count.
            
            reconstructed.topology =
                as.adj.matrix(reconstruction, models = m)[[m]]
            flag = TRUE;

            for (j in 1:nrow(reconstructed.topology)) {
                for (k in 1:ncol(reconstructed.topology)) {
                    if (reconstructed.topology[j, k] !=
                        curr.adj.matrix[j, k]) {
                        flag = FALSE
                        next()
                    }
                }
            }

            if (flag == TRUE) {
                overall.confidence[[m]] = overall.confidence[[m]] + 1
                overall.frequency[[m]] = overall.confidence[[m]] / nboot
            }
        }
    }
    
    ## Save the edge confidence and the frequency of the bootstrap
    ## adj.matrix.
    
    for (m in names(as.models(reconstruction))) {
        
        curr.adj.matrix = as.adj.matrix(reconstruction, models = m)[[m]];
        
        ## Save the edge confidence.
        
        curr.bootstrap.matrix =
            bootstrap.adj.matrix[[m]][-1, -1];
        curr.edge.confidence =
            array(0, c(ncol(curr.bootstrap.matrix), nrow(curr.bootstrap.matrix)))
        colnames(curr.edge.confidence) =
            colnames(curr.bootstrap.matrix);
        rownames(curr.edge.confidence) =
            rownames(curr.bootstrap.matrix);
        
        for (i in 1:ncol(curr.bootstrap.matrix)) {
            for (j in 1:nrow(curr.bootstrap.matrix)) {
                curr.edge.confidence[i, j] =
                    (curr.adj.matrix[i, j] *
                         as.numeric(curr.bootstrap.matrix[i, j])) / nboot
            }
        }
        bootstrap.edge.confidence[[m]] = curr.edge.confidence
        
        ## Save the frequency of the bootstrap adj.matrix.
        
        curr.bootstrap.matrix = bootstrap.adj.matrix[[m]];
        curr.adj.matrix.frequency =
            array(0, c(ncol(curr.bootstrap.matrix), nrow(curr.bootstrap.matrix)))
        colnames(curr.adj.matrix.frequency) =
            colnames(curr.bootstrap.matrix);
        rownames(curr.adj.matrix.frequency) =
            rownames(curr.bootstrap.matrix);
        
        for (i in 1:ncol(curr.bootstrap.matrix)) {
            for (j in 1:nrow(curr.bootstrap.matrix)) {
                curr.adj.matrix.frequency[i, j] =
                    as.numeric(as.numeric(curr.bootstrap.matrix[i, j])) / nboot
            }
        }
        bootstrap.adj.matrix.frequency[[m]] = curr.adj.matrix.frequency
    }
    
    ## Save the statistics of the bootstrap.
    
    for (m in names(as.models(reconstruction))) {
        if (command == "non-parametric") {
            
            bootstrap.statistics[[m]]$npb$bootstrap.results =
                bootstrap.results[[m]]
            
            bootstrap.statistics[[m]]$npb$bootstrap.adj.matrix =
                list(count = bootstrap.adj.matrix[[m]],
                     frequency = bootstrap.adj.matrix.frequency[[m]])
            
            bootstrap.statistics[[m]]$npb$bootstrap.edge.confidence =
                bootstrap.edge.confidence[[m]]
            
            bootstrap.statistics[[m]]$npb$overall.confidence =
                list(count = overall.confidence[[m]],
                     frequency = overall.frequency[[m]])
            
            bootstrap.statistics[[m]]$npb$bootstrap.settings =
                list(type = command, nboot = nboot)
        } else if (command == "parametric") {
            
            bootstrap.statistics[[m]]$pb$bootstrap.results =
                bootstrap.results[[m]]
            
            bootstrap.statistics[[m]]$pb$bootstrap.adj.matrix =
                list(count = bootstrap.adj.matrix[[m]],
                     frequency = bootstrap.adj.matrix.frequency[[m]])
            
            bootstrap.statistics[[m]]$pb$bootstrap.edge.confidence =
                bootstrap.edge.confidence[[m]]
            
            bootstrap.statistics[[m]]$pb$overall.confidence =
                list(count = overall.confidence[[m]],
                     frequency = overall.frequency[[m]])
            
            bootstrap.statistics[[m]]$pb$bootstrap.settings =
                list(type = command, nboot = nboot)
        } else if (command == "statistical") {
            
            bootstrap.statistics[[m]]$sb$bootstrap.results =
                bootstrap.results[[m]]
            
            bootstrap.statistics[[m]]$sb$bootstrap.adj.matrix =
                list(count = bootstrap.adj.matrix[[m]],
                     frequency = bootstrap.adj.matrix.frequency[[m]])
            
            bootstrap.statistics[[m]]$sb$bootstrap.edge.confidence =
                bootstrap.edge.confidence[[m]]
            
            bootstrap.statistics[[m]]$sb$overall.confidence =
                list(count = overall.confidence[[m]],
                     frequency = overall.frequency[[m]])
            
            bootstrap.statistics[[m]]$sb$bootstrap.settings =
                list(type = command, nboot = nboot)
        }
    }
    
    ## Save the execution time of the bootstrap.
    
    if (command == "non-parametric") {
        bootstrap.statistics$npb$execution.time = (proc.time() - ptm)
    } else if (command == "parametric") {
        bootstrap.statistics$pb$execution.time = (proc.time() - ptm)
    } else if (command == "statistical") {
        bootstrap.statistics$sb$execution.time = (proc.time() - ptm)
    }
    
    return(bootstrap.statistics)
}


#### end of file -- edmonds.bootstrap.R
