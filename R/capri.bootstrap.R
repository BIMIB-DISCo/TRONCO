#### bootstrap.capri.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#perform non-parametric or parametric bootstrap to evalutate the confidence of the reconstruction
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#hypotheses: a set of hypotheses referring to the dataset
#command.capri: type of search, either hill climbing (hc) or tabu (tabu)
#do.boot: should I perform bootstrap? Yes if TRUE, no otherwise
#nboot.capri: integer number (greater than 0) of bootstrap sampling to be performed
#pvalue: pvalue for the tests (value between 0 and 1)
#reconstructed.topology.pf: previously reconstructed prima facie topology (before any bootstrap)
#reconstructed.topology.bic: previously reconstructed topology (before any bootstrap)
#command: should I perform non-parametric or parametric bootstrap?
#estimated.marginal.probabilities.pf: estimated marginal probabilities of the events given the selected error rates for the prima facie topology
#estimated.conditional.probabilities.pf: estimated conditional probabilities of the events given the selected error rates for the prima facie topology
#parents.pos.pf: positions of the parents for prima facie
#error.rates.pf: selected error rates to be used if the bootstrap is "parametric" for the prima facie topology
#estimated.marginal.probabilities.bic: estimated marginal probabilities of the events given the selected error rates for the causal topology
#estimated.conditional.probabilities.bic: estimated conditional probabilities of the events given the selected error rates for the causal topology
#parents.pos.bic: positions of the parents for bic
#error.rates.bic: selected error rates to be used if the bootstrap is "parametric" for the causal topology
#nboot: number of bootstrap resampling to be performed
#RETURN:
#bootstrap.statistics: statistics of the bootstrap
bootstrap.capri <- function(dataset, 
                            hypotheses, 
                            command.capri, 
                            do.boot, nboot.capri, 
                            pvalue, reconstructed.topology.pf, 
                            reconstructed.topology.bic, 
                            command = "non-parametric", 
                            estimated.marginal.probabilities.pf, 
                            estimated.conditional.probabilities.pf, 
                            parents.pos.pf, 
                            error.rates.pf, 
                            estimated.marginal.probabilities.bic, 
                            estimated.conditional.probabilities.bic, 
                            parents.pos.bic, 
                            error.rates.bic, 
                            nboot, 
                            REGULARIZATION) 
{
    #structure to save the statistics of the bootstrap
    bootstrap.adj.matrix.pf = array(0, c(ncol(dataset) + 1, ncol(dataset) + 1))
    colnames(bootstrap.adj.matrix.pf) = c("None", colnames(dataset))
    rownames(bootstrap.adj.matrix.pf) = c("None", colnames(dataset))
    
    bootstrap.adj.matrix.bic = array(0, c(ncol(dataset) + 1, ncol(dataset) + 1))
    colnames(bootstrap.adj.matrix.bic) = c("None", colnames(dataset))
    rownames(bootstrap.adj.matrix.bic) = c("None", colnames(dataset))
    
    #structure to save the results of the bootstrap
    bootstrap.results.pf = array(list(-1), c(nboot, ncol(dataset)))
    colnames(bootstrap.results.pf) = colnames(dataset)

    bootstrap.results.bic = array(list(-1), c(nboot, ncol(dataset)))
    colnames(bootstrap.results.bic) = colnames(dataset)
    
    #perform nboot bootstrap resampling
  
  
    # create a progress bar
    flush.console()

    pb <- txtProgressBar(1, nboot, style = 3)
  
    for (num in 1:nboot) {
        
        #start the progress bar
        setTxtProgressBar(pb, num)
        
        #performed the bootstrapping procedure
        if(command == "non-parametric") {
            
            #perform the sampling for the current step of bootstrap
            samples = sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE)
            
            #perform the reconstruction on the bootstrapped dataset
            check.data = check.dataset(dataset[samples, ], FALSE)
            
            #if the reconstruction was performed without errors
            if(check.data$is.valid == TRUE) {
                bootstrapped.dataset = check.data$dataset
                curr.hypotheses = hypotheses

                if(!is.na(check.data$invalid.events$removed.events[1]) && !is.na(curr.hypotheses)) {
                    curr.hypotheses$num.hypotheses = hypotheses$num.hypotheses - 
                                                        length(which(unique(hypotheses$hlist[,"cause"]) %in% 
                                                        colnames(dataset)[check.data$invalid.events$removed.events]))
                }

                bootstrapped.hypotheses = curr.hypotheses
                bootstrapped.topology = capri.fit(bootstrapped.dataset, 
                                                  bootstrapped.hypotheses,
                                                  command.capri,
                                                  REGULARIZATION
                                                  do.boot,
                                                  nboot.capri,
                                                  pvalue,
                                                  3,
                                                  TRUE,
                                                  12345,
                                                  FALSE,
                                                  TRUE)
                
                #set the reconstructed causal edges
                parents.pos.pf = array(list(), c(ncol(bootstrapped.topology$data), 1))
                parents.pos.bic = array(list(), c(ncol(bootstrapped.topology$data),1))
                
                for(i in 1:ncol(bootstrapped.topology$data)) {
                    for(j in 1:ncol(bootstrapped.topology$data)) {
                        if(i != j && bootstrapped.topology$adj.matrix$adj.matrix.pf[i, j] == 1) {
                            parents.pos.pf[j, 1] = list(c(unlist(parents.pos.pf[j, 1]), i))
                        }
                        if(i != j && bootstrapped.topology$adj.matrix$adj.matrix.fit[i, j] == 1) {
                            parents.pos.bic[j, 1] = list(c(unlist(parents.pos.bic[j, 1]), i))
                        }
                    }
                }

                parents.pos.pf[unlist(lapply(parents.pos.pf, is.null))] = list(-1)
                parents.pos.bic[unlist(lapply(parents.pos.bic, is.null))] = list(-1)
                
                #get the matched edge in the reconstruction
                matched.idx.pf = match(colnames(bootstrapped.topology$data), colnames(bootstrap.results.pf))
                matched.idx.bic = match(colnames(bootstrapped.topology$data), colnames(bootstrap.results.bic))
                
                #if an event has no match, it means it has been merged and I discard it
                parents.pos.pf = parents.pos.pf[!is.na(matched.idx.pf)];
                parents.pos.bic = parents.pos.bic[!is.na(matched.idx.bic)];
                matched.idx.pf = matched.idx.pf[!is.na(matched.idx.pf)];
                matched.idx.bic = matched.idx.bic[!is.na(matched.idx.bic)];
                
                #save the results
                bootstrap.results.pf[num, matched.idx.pf] = parents.pos.pf
                bootstrap.results.bic[num, matched.idx.bic] = parents.pos.bic
            }

        } else if(command=="parametric") {

            if(num == 1) {

                #define the possible samples given the current number of events
                possible.strings = 2 ^ ncol(dataset)
                err = ""
                message = "Too many events! Parametric bootstrastap can not be performed."
                err <- tryCatch(curr.dataset = suppressWarnings(array(0, c(possible.strings, ncol(dataset)))),
                                                                error = function(e) err <- message)
                
                if(toString(err) == message) {
                    stop(err, call. = FALSE)
                }
                
                for (i in 1:possible.strings) {
                    curr.dataset[i, ] = decimal.to.binary.dag(i - 1, ncol(dataset))
                }

                colnames(curr.dataset) = colnames(dataset)

                #define the samples distribution induced by the topology
                samples.probabilities.pf = estimate.dag.samples(curr.dataset,
                                                                reconstructed.topology.pf,
                                                                estimated.marginal.probabilities.pf,
                                                                estimated.conditional.probabilities.pf,
                                                                parents.pos.pf,
                                                                error.rates.pf)

                samples.probabilities.bic = estimate.dag.samples(curr.dataset,
                                                                 reconstructed.topology.bic,
                                                                 estimated.marginal.probabilities.bic,
                                                                 estimated.conditional.probabilities.bic,
                                                                 parents.pos.bic,
                                                                 error.rates.bic)

            }

            #perform the sampling for the current step of bootstrap
            samples.pf <- suppressWarnings(sample(1:nrow(curr.dataset),
                                                  size=nrow(dataset),
                                                  replace=TRUE,
                                                  prob=samples.probabilities.pf))

            samples.bic <- suppressWarnings(sample(1:nrow(curr.dataset),
                                                   size=nrow(dataset),
                                                   replace=TRUE,
                                                   prob=samples.probabilities.bic))
            
            #perform the reconstruction on the bootstrapped dataset
            check.data.pf = check.dataset(curr.dataset[samples.pf, ], FALSE)
            check.data.bic = check.dataset(curr.dataset[samples.bic, ], FALSE)
            
            #if the reconstruction was performed without errors for the prima facie topology
            if(check.data.pf$is.valid == TRUE) {
                bootstrapped.dataset = check.data.pf$dataset
                curr.hypotheses = hypotheses

                if(!is.na(check.data.pf$invalid.events$removed.events[1]) && !is.na(curr.hypotheses)) {
                    curr.hypotheses$num.hypotheses = hypotheses$num.hypotheses - 
                                                     length(which(unique(hypotheses$hlist[,"cause"]) %in% 
                                                     colnames(dataset)[check.data$invalid.events$removed.events]))
                }
                
                bootstrapped.hypotheses = curr.hypotheses
                bootstrapped.topology = capri.fit(bootstrapped.dataset, 
                                                  bootstrapped.hypotheses,
                                                  command.capri,
                                                  REGULARIZATION
                                                  do.boot,
                                                  nboot.capri,
                                                  pvalue,
                                                  3,
                                                  TRUE,
                                                  12345,
                                                  FALSE,
                                                  TRUE)
                
                #set the reconstructed causal edges
                parents.pos.pf = array(list(), c(ncol(bootstrapped.topology$data), 1))

                for(i in 1:ncol(bootstrapped.topology$data)) {
                    for(j in 1:ncol(bootstrapped.topology$data)) {
                        if(i != j && bootstrapped.topology$adj.matrix$adj.matrix.pf[i, j] == 1) {
                            parents.pos.pf[j, 1] = list(c(unlist(parents.pos.pf[j, 1]), i))
                        }
                    }
                }

                parents.pos.pf[unlist(lapply(parents.pos.pf, is.null))] = list(-1)

                #get the matched edge in the reconstruction
                matched.idx.pf = match(colnames(bootstrapped.topology$data), colnames(bootstrap.results.pf))

                #if an event has no match, it means it has been merged and I discard it
                parents.pos.pf = parents.pos.pf[!is.na(matched.idx.pf)]
                matched.idx.pf = matched.idx.pf[!is.na(matched.idx.pf)]

                #save the results
                bootstrap.results.pf[num, matched.idx.pf] = parents.pos.pf;
            }

            #if the reconstruction was performed without errors for the causal topology
            if(check.data.bic$is.valid == TRUE) {
                bootstrapped.dataset = check.data.bic$dataset
                curr.hypotheses = hypotheses

                if(!is.na(check.data.bic$invalid.events$removed.events[1]) && !is.na(curr.hypotheses)) {
                    curr.hypotheses$num.hypotheses = hypotheses$num.hypotheses - 
                                                     length(which(unique(hypotheses$hlist[,"cause"]) %in%
                                                      colnames(dataset)[check.data$invalid.events$removed.events]))
                }

                bootstrapped.hypotheses = curr.hypotheses
                bootstrapped.topology = capri.fit(bootstrapped.dataset, 
                                                  bootstrapped.hypotheses,
                                                  command.capri,
                                                  REGULARIZATION
                                                  do.boot,
                                                  nboot.capri,
                                                  pvalue,
                                                  3,
                                                  TRUE,
                                                  12345,
                                                  FALSE,
                                                  TRUE)

                #set the reconstructed causal edges
                parents.pos.bic = array(list(), c(ncol(bootstrapped.topology$data), 1))

                for(i in 1:ncol(bootstrapped.topology$data)) {
                    for(j in 1:ncol(bootstrapped.topology$data)) {
                        if(i != j && bootstrapped.topology$adj.matrix$adj.matrix.fit[i, j] == 1) {
                            parents.pos.bic[j, 1] = list(c(unlist(parents.pos.bic[j, 1]), i))
                        }
                    }
                }

                parents.pos.bic[unlist(lapply(parents.pos.bic, is.null))] = list(-1)

                #get the matched edge in the reconstruction
                matched.idx.bic = match(colnames(bootstrapped.topology$data), colnames(bootstrap.results.bic))

                #if an event has no match, it means it has been merged and I discard it
                parents.pos.bic = parents.pos.bic[!is.na(matched.idx.bic)]
                matched.idx.bic = matched.idx.bic[!is.na(matched.idx.bic)]

                #save the results
                bootstrap.results.bic[num, matched.idx.bic] = parents.pos.bic
            }
        }
    }
  
    # close progress bar
    close(pb)
  
    #set the statistics of the bootstrap
    for(i in 2:ncol(bootstrap.adj.matrix.pf)) {

        curr.result.pf = bootstrap.results.pf[, i - 1]
        curr.result.bic = bootstrap.results.bic[, i - 1]

        for(j in 1:length(curr.result.pf)) {
            curr.pf = curr.result.pf[[j]]
            
            for(k in 1:length(curr.pf)) {
                if(curr.pf[k] == -1) {
                    bootstrap.adj.matrix.pf[1, i] = bootstrap.adj.matrix.pf[1, i] + 1
                } else {
                    bootstrap.adj.matrix.pf[curr.pf[k] + 1, i] = bootstrap.adj.matrix.pf[curr.pf[k] + 1, i] + 1
                }
            }
            curr.bic = curr.result.bic[[j]]

            for(k in 1:length(curr.bic)) {
                if(curr.bic[k] == -1) {
                    bootstrap.adj.matrix.bic[1, i] = bootstrap.adj.matrix.bic[1, i] + 1
                } else {
                    bootstrap.adj.matrix.bic[curr.bic[k] + 1, i] = bootstrap.adj.matrix.bic[curr.bic[k] + 1, i] + 1
                }
            }
        }
    }

    #evalutate the overall confidence
    overall.confidence.pf = 0
    overall.confidence.bic = 0
    cont = 0

    for(i in 1:nrow(bootstrap.results.pf)) {
        curr.adj.matrix.pf = array(0, c(ncol(dataset), ncol(dataset)))
        curr.adj.matrix.bic = array(0, c(ncol(dataset), ncol(dataset)))
        
        for(j in 1:ncol(bootstrap.results.pf)) {
            curr.result.pf = bootstrap.results.pf[i, j]
            curr.result.bic = bootstrap.results.bic[i, j]

            for(k in 1:length(curr.result.pf)) {
                curr.pf = curr.result.pf[[k]]
                
                for(l in 1:length(curr.pf)) {
                    if(curr.pf[l] != -1) {
                        curr.adj.matrix.pf[curr.pf[l], j] = 1
                    }
                }
                curr.bic = curr.result.bic[[k]]
                for(l in 1:length(curr.bic)) {
                    if(curr.bic[l] != -1) {
                        curr.adj.matrix.bic[curr.bic[l],j] = 1
                    }
                }
            }
        }
        if(sum(reconstructed.topology.pf - curr.adj.matrix.pf) == 0) {
            overall.confidence.pf = overall.confidence.pf + 1
        }
        if(sum(reconstructed.topology.bic - curr.adj.matrix.bic) == 0) {
            overall.confidence.bic = overall.confidence.bic + 1
        }
    }

    #save the reconstructed topologies
    reconstructed.topology = list(reconstructed.topology.pf = reconstructed.topology.pf,
                                  reconstructed.topology.bic = reconstructed.topology.bic)
    
    #save the edge confidence
    edge.confidence.pf = (reconstructed.topology.pf * bootstrap.adj.matrix.pf[-1, -1]) / nboot
    edge.confidence.bic = (reconstructed.topology.bic * bootstrap.adj.matrix.bic[-1, -1])/nboot
    edge.confidence = list(edge.confidence.pf=edge.confidence.pf, edge.confidence.bic=edge.confidence.bic)
    
    #save the confidence from the bootstrap
    confidence.pf = list(overall.value.pf = overall.confidence.pf,
                         overall.frequency.pf = overall.confidence.pf / nboot,
                         bootstrap.values.pf = bootstrap.adj.matrix.pf[, -1],
                         bootstrap.frequencies.pf = bootstrap.adj.matrix.pf[, -1] / nboot)

    confidence.bic = list(overall.value.bic = overall.confidence.bic,
                          overall.frequency.bic = overall.confidence.bic / nboot,
                          bootstrap.values.bic = bootstrap.adj.matrix.bic[, -1],
                          bootstrap.frequencies.bic = bootstrap.adj.matrix.bic[, -1] / nboot)
    
    confidence = list(confidence.pf = confidence.pf, confidence.bic = confidence.bic)   
    
    #save the settings of the bootstrap
    bootstrap.settings = list(type = command, nboot = nboot)
    
    #structure to save the results
    bootstrap.statistics = list(reconstructed.topology = reconstructed.topology,
                                confidence = confidence,
                                edge.confidence = edge.confidence,
                                bootstrap.settings = bootstrap.settings)

    return(bootstrap.statistics)
}

#### end of file -- bootstrap.capri.R


#### decimal.to.binary.dag.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#convert an integer decimal number to binary
#INPUT:
#num.decimal: decimal integer to be converted
#num.bits: number of bits to be used
#RETURN:
#num.binary: binary conversion of num.decimal
decimal.to.binary.dag <- function(num.decimal, num.bits) {
    
    #structure where to save the result
    num.binary = rep(0, num.bits)
    
    #convert the integer decimal number to binary
    pos = 0
    while(num.decimal > 0) {
        
        #compute the value of the current step
        num.binary[num.bits-pos] = num.decimal %% 2;
        
        #divide the number by 2 for the next iteration
        num.decimal = num.decimal %/% 2
        pos = pos + 1
    }
    return(num.binary)
}

#### end of file -- decimal.to.binary.dag.R


#### estimate.dag.samples.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the probability of observing each sample in the dataset given the reconstructed topology
#INPUT:
#dataset: a valid dataset
#reconstructed.topology: the reconstructed topology
#estimated.marginal.probabilities: estimated marginal probabilities of the events
#estimated.conditional.probabilities: estimated conditional probabilities of the events
#parents.pos: position of the parents of each node
#error.rates: error rates for false positives and false negatives
#RETURN:
#probabilities: probability of each sample
estimate.dag.samples = function(dataset, 
                                reconstructed.topology, 
                                estimated.marginal.probabilities, 
                                estimated.conditional.probabilities, 
                                parents.pos, error.rates) 
{
    #structure where to save the probabilities of the samples
    probabilities = array(-1, c(nrow(dataset), 1))
    
    #compute the position of the latest parent and its conditional probability for each node
    last.parent.pos = array(-1, c(nrow(parents.pos), 1))
    curr.estimated.conditional.probabilities = array(1, c(nrow(estimated.conditional.probabilities), 1))
    
    for (i in 1:length(parents.pos)) {
        if(length(parents.pos[[i, 1]]) != 1 || parents.pos[[i, 1]] != -1) {
            curr.last.parent = which.min(estimated.marginal.probabilities[parents.pos[[i, 1]], 1])
            last.parent.pos[i, 1] = parents.pos[[i, 1]][curr.last.parent[1]]
            curr.estimated.conditional.probabilities[i, 1] = estimated.conditional.probabilities[[i, 1]][curr.last.parent[1]]
        }
    }
    
    #topological properties:
    #1. progression number
    #2. latest parent
    #3. level in the progression
    
    topology.structure = array(0, c(nrow(reconstructed.topology), 3))
    
    #go through the subtrees within the topology
    progression.count = 0
    
    for (i in 1:nrow(reconstructed.topology)) {
        
        #if node i has no parents, it is a root
        if(length(which(reconstructed.topology[, i] == 1)) == 0) {
            progression.count = progression.count + 1
            level = 1
            
            #set the parameters for the root
            topology.structure[i,1] = progression.count
            topology.structure[i,2] = -1
            topology.structure[i,3] = level
            curr.node = i
            
            #go through this progression
            while (length(curr.node) > 0) {
                
                #move to the next level
                level = level + 1
                new.node = vector()
                
                for (j in 1:length(curr.node)) {
                    curr.new.node = which(reconstructed.topology[curr.node[j], ] == 1)
                    
                    if(length(curr.new.node) > 0) {
                        new.node = c(new.node,curr.new.node)
                        
                        for (k in 1:length(curr.new.node)) {
                            
                            #number of the current subprogression
                            topology.structure[curr.new.node[k], 1] = progression.count
                            
                            #parent of the current node
                            if(last.parent.pos[curr.new.node[k], 1] == curr.node[j]) {
                                topology.structure[curr.new.node[k], 2] = curr.node[j]
                            }
                            
                            #level of this node
                            topology.structure[curr.new.node[k], 3] = level
                        }
                    }
                }
                curr.node = new.node
            }
        }
    }

    #go through the dataset and evalutate the probability of each sample
    for (i in 1:nrow(dataset)) {
        sample.probability = 1

        for (j in 1:progression.count) {
            
            #probability of this subprogression (without any knowledge, I set it to 1)
            curr.sample.probability = 1
            
            #entries referring to this subprogression
            curr.entry = which(topology.structure[, 1] == j)
            
            #samples of each element of this subprogression
            curr.sample = dataset[i, curr.entry]
            
            #parents of each element of this subprogression
            curr.parents = topology.structure[curr.entry, 2]
            
            #level of each element of this subprogression
            curr.levels = topology.structure[curr.entry, 3]
            
            #set the probability as the one of the root of this progression
            curr.sample.probability = curr.sample.probability * estimated.marginal.probabilities[curr.entry[which(curr.levels == 1, arr.ind = TRUE)], 1]
            
            #set the maximum level of this subprogression
            max.level = curr.levels[which.max(curr.levels)]

            #if I have at least one event in this sample
            if(length(curr.sample[curr.sample == 1]) > 0) {

                #visit the nodes starting from the lower level
                is.valid = TRUE

                for (k in max.level:1) {
                    curr.level.nodes = which(curr.levels == k, arr.ind=TRUE)

                    #if I'm not on a root
                    if(k > 1) {
                        curr.level.samples = curr.sample[curr.level.nodes]

                        #if I have at least one event at this level
                        if(length(curr.level.samples[curr.level.samples == 1]) > 0) {
                            
                            #I can not have a child without its parent
                            curr.level.parent = curr.parents[curr.level.nodes]
                            
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i, curr.level.parent[p]] == 0 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    is.valid = FALSE
                                    break
                                }
                            }
                        }

                        #if the sample is valid
                        if(is.valid == TRUE) {
                            
                            #add the probability of each edge
                            curr.level.parent = curr.parents[curr.level.nodes]
                            
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i, curr.level.parent[p]] == 1 && dataset[i,curr.entry[curr.level.nodes[p]]] == 0) {
                                    curr.sample.probability = curr.sample.probability * (1 - curr.estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1])
                                } else if(dataset[i, curr.level.parent[p]]==1 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    curr.sample.probability = curr.sample.probability * curr.estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1]
                                }
                            }
                        }
                    }

                    if(is.valid == FALSE) {
                        curr.sample.probability = 0
                        break
                    }
                }
                if(is.valid == FALSE) {
                    sample.probability = 0
                    break
                }
            } else {
                # ..if this sample has no events for this progression
                curr.sample.probability = 1 - curr.sample.probability
            }

            #update the probability of the topology with the one of this sample
            sample.probability = sample.probability * curr.sample.probability

            if(sample.probability == 0) {
                break
            }
        }
        probabilities[i, 1] = sample.probability
    }

    #correct the estimation by the error rates
    errors.matrix = array(0, c(nrow(probabilities), nrow(dataset)))
    for (i in 1:nrow(probabilities)) {
        for (j in 1:nrow(dataset)) {
            curr.sample.x = as.numeric(dataset[i, ])
            curr.sample.y = as.numeric(dataset[j, ])
            errors.matrix[i, j] = (1 - error.rates$error.fp) ^ ((1 - curr.sample.x) %*% (1 - curr.sample.y)) *
                                error.rates$error.fp ^ ((1 - curr.sample.x) %*% curr.sample.y) *
                                (1 - error.rates$error.fn) ^ (curr.sample.x %*% curr.sample.y) *
                                error.rates$error.fn ^ (curr.sample.x %*% (1 - curr.sample.y))
        }
    }
    probabilities[, 1] = as.numeric(as.vector(probabilities) %*% errors.matrix)
    return(probabilities)
}

#### end of file -- estimate.dag.samples.R


