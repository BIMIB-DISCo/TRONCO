#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# perform non-parametric or statistical bootstrap to evalutate the confidence of the reconstruction
# @title bootstrap.capri
# @param reconstruction Result of a previous reconstruction
# @param command should I perform non-parametric or statistical bootstrap?
# @param nboot.algorithm integer number (greater than 0) of bootstrap sampling to be performed
# @param cores.ratio Percentage of cores to use
# @return bootstrap.statistics: statistics of the bootstrap
#
bootstrap <- function(reconstruction, 
                      command = "non-parametric",
                      nboot = 100,
                      cores.ratio = 1,
                      silent = FALSE) {
    

    parameters = as.parameters(reconstruction)
    type = parameters$algorithm

    dataset = as.genotypes(reconstruction)

    bootstrap.statistics = list()
    if (!is.null(reconstruction$bootstrap)) {
        bootstrap.statistics = reconstruction$bootstrap
    }
    
    ## Get CAPRI parameters.
    
    if (type == 'CAPRI') {
        command.capri = parameters$command
        hypotheses = NA
        if (nhypotheses(reconstruction) > 0) {
            hypotheses = reconstruction$hypotheses
        }
    }

    ## Get CAPRESE parameters.

    if (type == 'CAPRESE') {
        lambda = parameters$lambda
    }

    ## Get other parameters.

    if (type %in% c('CAPRI', 'EDMONDS', 'CHOW_LIU', 'PRIM')) {
        regularization = parameters$regularization
        do.boot = parameters$do.boot
        nboot.algorithm = parameters$nboot
        pvalue = parameters$pvalue
        min.boot = parameters$min.boot
        min.stat = parameters$min.stat
        boot.seed = parameters$boot.seed
    }

    ## Start the clock to measure the execution time

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
    if (!silent) {
        cat("\tExpected completion in approx.",
            format(.POSIXct(expected.execution.time, tz="GMT"),
                   "%Hh:%Mm:%Ss"),
            "\n")
    }


    cl = makeCluster(cores)    
    registerDoParallel(cl)
    if (!silent) {
        cat('\tUsing', cores, 'cores via "parallel" \n')
    }
    ## Perform nboot bootstrap resampling
    
    r = foreach(num = 1:nboot ) %dopar% {    
            
        ## Performed the bootstrapping procedure.
        
        curr.reconstruction = list()

        if (command == "non-parametric") {
            
            ## Perform the sampling for the current step of bootstrap.
            
            samples = sample(1:nrow(dataset), 
                             size = nrow(dataset),
                             replace = TRUE)
            curr.reconstruction$genotypes = dataset[samples,]


        } else if (command == "statistical") {
            curr.reconstruction$genotypes = reconstruction$genotypes;
        }

        curr.reconstruction$annotations = reconstruction$annotations
        curr.reconstruction$types = reconstruction$types

        if (type == 'CAPRI') {
            curr.reconstruction$hypotheses = hypotheses
            bootstrapped.topology =
                tronco.capri(curr.reconstruction,
                             command.capri, 
                             regularization,
                             do.boot,
                             nboot.algorithm,
                             pvalue,
                             min.boot,
                             min.stat,
                             boot.seed,
                             silent = TRUE)
        } else if (type == 'CAPRESE') {
            bootstrapped.topology =
                tronco.caprese(curr.reconstruction,
                               lambda, 
                               silent = TRUE)
        } else if (type == 'CHOW_LIU') {
            bootstrapped.topology = 
                tronco.mst.chowliu(curr.reconstruction,
                                   regularization,
                                   do.boot,
                                   nboot.algorithm,
                                   pvalue,
                                   min.boot,
                                   min.stat,
                                   boot.seed,
                                   silent = TRUE)
        } else if (type == 'PRIM') {
            bootstrapped.topology =
                tronco.mst.prim(curr.reconstruction,
                                regularization,
                                do.boot,
                                nboot.algorithm,
                                pvalue,
                                min.boot,
                                min.stat,
                                boot.seed,
                                silent = TRUE)
        } else if (type == 'EDMONDS') {
            bootstrapped.topology =
                tronco.mst.edmonds(curr.reconstruction,
                                   regularization,
                                   do.boot,
                                   nboot.algorithm,
                                   pvalue,
                                   min.boot,
                                   min.stat,
                                   boot.seed,
                                   silent = TRUE)
        }
            
        curr.reconstruction = bootstrapped.topology
        
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
    if (!silent) {
        cat("\tReducing results\n")
    }
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


# convert an integer decimal number to binary
# @title decimal.to.binary.dag
# @param num.decimal decimal integer to be converted
# @param num.bits number of bits to be used
# @return num.binary: binary conversion of num.decimal
#
decimal.to.binary.dag <- function(num.decimal, num.bits) {
    
    ## Structure where to save the result.
    
    num.binary = rep(0, num.bits)
    
    ## Convert the integer decimal number to binary.
    
    pos = 0
    while (num.decimal > 0) {
        
        ## Compute the value of the current step.
        
        num.binary[num.bits-pos] = num.decimal %% 2;
        
        ## Divide the number by 2 for the next iteration.
        
        num.decimal = num.decimal %/% 2
        pos = pos + 1
    }
    return(num.binary)
}


#### end of file -- capri.bootstrap.R
