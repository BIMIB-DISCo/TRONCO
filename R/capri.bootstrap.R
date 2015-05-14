#### bootstrap.capri.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# perform non-parametric or parametric bootstrap to evalutate the confidence of the reconstruction
# INPUT:
# dataset: a dataset describing a progressive phenomenon
# hypotheses: a set of hypotheses referring to the dataset
# command.capri: type of search, either hill climbing (hc) or tabu (tabu)
# regularization: regularizators to be used for the likelihood fit
# do.boot: should I perform bootstrap? Yes if TRUE, no otherwise
# nboot.capri: integer number (greater than 0) of bootstrap sampling to be performed
# pvalue: pvalue for the tests (value between 0 and 1)
# min.boot: minimum number of bootstrapping to be performed
# min.stat: should I keep bootstrapping untill I have nboot valid values?
# boot.seed: seed to be used for the sampling
# do.estimation: should I perform the estimation of the error rates and probabilities?
# silent: should I be verbose?
# command: should I perform non-parametric or parametric bootstrap?
# nboot: number of bootstrap resampling to be performed
# RETURN:
# bootstrap.statistics: statistics of the bootstrap
bootstrap.capri <- function(dataset, 
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
                            command = "non-parametric",
                            nboot) 
{
    
    # structure to save the results of the bootstrap
    curr.bootstrap.results = array(list(-1), c(nboot,nevents(reconstruct)))
    colnames(curr.bootstrap.results) = as.events(reconstruct)
    bootstrap.results = list()
    bootstrap.results[names(as.models(reconstruction))] = curr.bootstrap.results
    
    curr.bootstrap.adj.matrix = array(0, c(nevents(reconstruct)+1,nevents(reconstruct)+1))
    colnames(curr.bootstrap.adj.matrix) = c("None",as.events(reconstruct))
    rownames(curr.bootstrap.adj.matrix) = c("None",as.events(reconstruct))
    bootstrap.adj.matrix = list()
    bootstrap.adj.matrix[names(as.models(reconstruction))] = curr.bootstrap.adj.matrix
    
    curr.bootstrap.adj.matrix = array(0, c(nevents(reconstruct)+1,nevents(reconstruct)+1))
    colnames(curr.bootstrap.adj.matrix) = c("None",as.events(reconstruct))
    rownames(curr.bootstrap.adj.matrix) = c("None",as.events(reconstruct))
    bootstrap.adj.matrix.frequency = list()
    bootstrap.adj.matrix.frequency[names(as.models(reconstruction))] = curr.bootstrap.adj.matrix
    
    curr.edge.confidence = array(0, c(nevents(reconstruct),nevents(reconstruct)))
    colnames(curr.edge.confidence) = as.events(reconstruct)
    rownames(curr.edge.confidence) = as.events(reconstruct)
    bootstrap.edge.confidence = list()
    bootstrap.edge.confidence[names(as.models(reconstruction))] = curr.edge.confidence
    
    overall.confidence = list();
    overall.confidence[names(as.models(reconstruction))] = 0;
    overall.frequency = list();
    overall.frequency[names(as.models(reconstruction))] = 0;
    
    # create a progress bar
    flush.console()
    pb <- txtProgressBar(1, nboot, style = 3)
  
  	# perform nboot bootstrap resampling
    for (num in 1:nboot) {
        
        # update the progress bar
        setTxtProgressBar(pb, num)
        
        # performed the bootstrapping procedure
        if(command == "non-parametric") {
            
            # perform the sampling for the current step of bootstrap
            samples = sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE)
            bootstrapped.dataset = dataset[samples,]
            
            # perform the reconstruction on the bootstrapped dataset
            bootstrapped.topology = capri.fit(bootstrapped.dataset, 
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
                            				  silent)
            
            curr.reconstruction = bootstrapped.topology;

        }
        else if(command=="parametric") {
        		
        		if(num == 1) {
        			
        			# structure to save the samples probabilities
        			samples.probabilities = list();
        			
                # define the possible samples given the current number of events
                possible.strings = 2 ^ ncol(dataset)
                err = ""
                message = "Too many events in the dataset! Parametric bootstrastap can not be performed."
                err <- tryCatch(curr.dataset = suppressWarnings(array(0, c(possible.strings, ncol(dataset)))),
                                                                error = function(e) err <- message)
                
                if(toString(err) == message) {
                		stop(err, call. = FALSE)
                }
                
                for (i in 1:possible.strings) {
                		curr.dataset[i, ] = decimal.to.binary.dag(i - 1, ncol(dataset))
                }

                colnames(curr.dataset) = colnames(dataset)
                
                for (m in names(as.models(reconstruction))) {
                		samples.probabilities[m] = estimate.dag.samples(curr.dataset,
                                                                as.adj.matrix(reconstruction,model=m),
                                                                as.marginal.probs(reconstruction,model=m,type="fit"),
                                                                as.conditional.probs(reconstruction,model=m,type="fit"),
                                                                as.parents.pos(reconstruction,model=m),
                                                                as.error.rates(reconstruction,model=m))
                }

            }
            
            # perform the reconstruction for each model
            curr.reconstruction = reconstuction;
            curr.reconstruction$models = list();
            for (m in names(as.models(reconstruction))) {
            	
            	# perform the sampling for the current step of bootstrap and regularizator
            	samples = sample(1:nrow(curr.dataset), size = nrow(dataset), replace = TRUE, prob = samples.probabilities[[m]])
            	bootstrapped.dataset = dataset[samples,]
            
            	# perform the reconstruction on the bootstrapped dataset
            	bootstrapped.topology = capri.fit(bootstrapped.dataset, 
                            				  	  hypotheses,
                            				  	  command.capri, 
                            				  	  m,
                            				  	  do.boot,
                            				  	  nboot.capri,
                            				  	  pvalue,
                            				  	  min.boot,
                            				  	  min.stat,
                            				  	  boot.seed,
                            				  	  do.estimation,
                            				      silent)
                            				      
                # save the results for this model
                curr.reconstruction$models[m] = as.models(bootstrapped.topology,models=m)
            	
            }
            
        }
        
        # set the reconstructed selective advantage edges
        for (m in names(as.models(curr.reconstruction))) {
        	
        	# get the parents pos
        	parents.pos = array(list(), c(nevents(curr.reconstruction), 1))
        	
        	curr.adj.matrix = as.adj.matrix(curr.reconstruction,model=m)
        	for(i in 1:nevents(curr.reconstruction)) {
        		for(j in 1:nevents(curr.reconstruction)) {
        			if(i!=j && curr.adj.matrix[i,j]==1) {
        				parents.pos[j, 1] = list(c(unlist(parents.pos[j,1]),i))
        			}
        		}
        	}
        	
        	parents.pos[unlist(lapply(parents.pos,is.null))] = list(-1)
        	
        	# get the matched edge in the reconstruction
        	matched.idx = match(as.events(curr.reconstruction),as.events(curr.reconstruction))
        	parents.pos = parents.pos[!is.na(matched.idx)]
        	matched.idx = matched.idx[!is.na(matched.idx)]
        	
        	# save the results
        	bootstrap.results[m][num, matched.idx] = parents.pos;
        	
        }
        
    }
  
    # close progress bar
    close(pb)
    
    # set the statistics of the bootstrap
    for (m in names(as.models(curr.reconstruction))) {
    	
    	curr.bootstrap.adj.matrix = bootstrap.adj.matrix[m]
    	
    	for(i in 2:ncol(curr.bootstrap.adj.matrix)) {
    		curr.result = bootstrap.results[m][,i-1]
    		
    		for(j in 1:length(curr.result)) {
    			
    			curr.val = curr.result[[j]]
    			
    			for(k in 1:length(curr.val)) {
    				if(curr.val[k] == -1) {
    					curr.bootstrap.adj.matrix[1,i] = curr.bootstrap.adj.matrix[1,i] + 1
    				}
    				else {
    					curr.bootstrap.adj.matrix[curr.val[k] + 1, i] = curr.bootstrap.adj.matrix[curr.val[k] + 1, i] + 1
    				}
    			}
    		
    		}
    		
    	}
    	
    	}
    
    # evalutate the overall confidence
    for (m in names(as.models(curr.reconstruction))) {
    	
    	curr.bootstrap.results = bootstrap.results[m]
    	
    	for(i in 1:nrow(curr.bootstrap.results)) {
    		
    		curr.adj.matrix = array(0, c(nevents(reconstruction),nevents(reconstruction)))
    		
    		for(j in 1:ncol(curr.bootstrap.results)) {
    			
    			curr.result = curr.bootstrap.results[i, j]
    			
    			for(k in 1:length(curr.result)) {
    				
    				curr.val = curr.result[[k]]
                
                	for(l in 1:length(curr.val)) {
                    	if(curr.val[l] != -1) {
                        	curr.adj.matrix[curr.val[l], j] = 1
                    	}
                	}
                	
    			}
    			
    		}
        	
        	# if I have a perfect match between the reconstructed topologies, increase the count
        	if(sum(as.adj.matrix(reconstruction,model=m) - curr.adj.matrix) == 0) {
            	overall.confidence[k] = overall.confidence.pf + 1
            	overall.frequency[k] = overall.confidence[k] / nboot
        	}
    		
    	}
    	
    }
    
    # save the edge confidence
    for (m in names(as.models(curr.reconstruction))) {
    	
    	curr.edge.confidence = (as.adj.matrix(curr.reconstruction,model=m) * bootstrap.adj.matrix[m][-1,-1]) / nboot
    	bootstrap.confidence[m] = curr.edge.confidence;
    	
    	curr.bootstrap.adj.matrix = bootstrap.adj.matrix[m][, -1] / nboot
    	bootstrap.adj.matrix.frequency[m] = curr.bootstrap.adj.matrix
    
    }
        
    # save the statistics of the bootstrap
    bootstrap.statistics = list(bootstrap.results = bootstrap.results,
    							bootstrap.adj.matrix = list(count = bootstrap.adj.matrix, frequency = bootstrap.adj.matrix.frequency),
    							bootstrap.edge.confidence = bootstrap.edge.confidence,
    							overall.confidence = list(count = overall.confidence, frequency = overall.frequency),
                                bootstrap.settings = list(type = command, nboot = nboot))

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


