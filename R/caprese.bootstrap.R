#### bootstrap.caprese.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#perform non-parametric or parametric bootstrap to evalutate the confidence of the reconstruction
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#lambda: shrinkage parameter (value in [0,1])
#reconstructed.topology: previously reconstructed topology (before any bootstrap)
#command: should I perform non-parametric or parametric bootstrap?
#estimated.marginal.probabilities: estimated marginal probabilities of the events given the selected error rates
#estimated.conditional.probabilities: estimated conditional probabilities of the events given the selected error rates
#error.rates: selected error rates to be used if the bootstrap is "parametric"
#nboot: number of bootstrap resampling to be performed
#RETURN:
#bootstrap.statistics: statistics of the bootstrap
"bootstrap.caprese" <- function(dataset, 
                            lambda,
                            do.estimation,
                            silent,
                            reconstruction, 
                            command = "non-parametric",
                            nboot = 100) 
{
	
	# structure to save the results of the bootstrap
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
    
    overall.confidence = list();
    overall.confidence[names(as.models(reconstruction))] = list(0);
    overall.frequency = list();
    overall.frequency[names(as.models(reconstruction))] = list(0);
    
    # reset the seed
    set.seed(NULL)
    
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
            
            curr.reconstruction = list()
            curr.reconstruction$genotypes = bootstrapped.dataset;
            curr.reconstruction$annotations = reconstruction$annotations;
            curr.reconstruction$types = reconstruction$types;
            curr.reconstruction$hypotheses = reconstruction$hypotheses;
            
            # perform the reconstruction on the bootstrapped dataset
            bootstrapped.topology = tronco.caprese(curr.reconstruction,
                            				  lambda, 
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
                		
                		# estimate the samples probabilities for each model
                		samples.probabilities[m] = estimate.dag.samples(curr.dataset,
                                                                as.adj.matrix(reconstruction,model=m)[[m]],
                                                                as.marginal.probs(reconstruction,model=m,type="fit")[[m]],
                                                                as.conditional.probs(reconstruction,model=m,type="fit")[[m]],
                                                                as.parents.pos(reconstruction,model=m)[[m]],
                                                                as.error.rates(reconstruction,model=m)[[m]])
                
                }

            }
            
            # perform the reconstruction for each model
            new.reconstruction = reconstuction;
            new.reconstruction$model = list();
            for (m in names(as.models(reconstruction))) {
            	
	            # perform the sampling for the current step of bootstrap and regularizator
	            samples = sample(1:nrow(curr.dataset), size = nrow(dataset), replace = TRUE, prob = samples.probabilities[[m]])
	            bootstrapped.dataset = dataset[samples,]
            
	            curr.reconstruction = list()
	            curr.reconstruction$genotypes = bootstrapped.dataset;
	            curr.reconstruction$annotations = reconstruction$annotations;
	            curr.reconstruction$types = reconstruction$types;
	            curr.reconstruction$hypotheses = reconstruction$hypotheses;
	            
            	# perform the reconstruction on the bootstrapped dataset
            	bootstrapped.topology = tronco.caprese(curr.reconstruction,
                            				  lambda, 
                            				  do.estimation,
                            				  silent)
                            				      
                # save the results for this model
                new.reconstruction$model[m] = as.models(bootstrapped.topology,models=m)
            	
            }
            curr.reconstruction = new.reconstruction;
            
        }
        
        # set the reconstructed selective advantage edges
        for (m in names(as.models(curr.reconstruction))) {
        	
	        	# get the parents pos
	        	parents.pos = array(list(), c(nevents(curr.reconstruction), 1))
	        	
	        	
	        	curr.adj.matrix = as.adj.matrix(curr.reconstruction,model=m)[[m]]
	        	for(i in 1:nevents(curr.reconstruction)) {
	        		for(j in 1:nevents(curr.reconstruction)) {
	        			if(i!=j && curr.adj.matrix[i,j]==1) {
	        				parents.pos[j, 1] = list(c(unlist(parents.pos[j,1]),i))
	        			}
	        		}
	        	}
	        	
	        	parents.pos[unlist(lapply(parents.pos,is.null))] = list(-1)
	        	
	        	# save the results
	        	bootstrap.results[[m]][num,] = parents.pos;
        	
        }
        
    }
  
    # close progress bar
    close(pb)
    
    # set the statistics of the bootstrap
    for (m in names(as.models(reconstruction))) {
    	
	    	curr.bootstrap.adj.matrix = bootstrap.adj.matrix[[m]]
	    	
	    	for(i in 2:ncol(curr.bootstrap.adj.matrix)) {
	    		
	    		curr.result = bootstrap.results[[m]][,i-1]
	    		
	    		for(j in 1:length(curr.result)) {
	    			
	    			curr.val = curr.result[[j]]
	    			
	    			for(k in 1:length(curr.val)) {
	    				if(length(curr.val[k])==1 && curr.val[k] == -1) {
	    					curr.bootstrap.adj.matrix[[1,i]] = curr.bootstrap.adj.matrix[[1,i]] + 1
	    				}
	    				else {
	    					curr.bootstrap.adj.matrix[[curr.val[k] + 1, i]] = curr.bootstrap.adj.matrix[[curr.val[k] + 1, i]] + 1
	    				}
	    			}
	    		
	    		}
	    		
	    	}
	    	
	    	bootstrap.adj.matrix[[m]] = curr.bootstrap.adj.matrix;
	    rownames(bootstrap.results[[m]]) =  paste("Iteration ",1:nrow(bootstrap.results[[m]]),sep="")
    	
    	}
    
    # evalutate the overall confidence
    for (m in names(as.models(reconstruction))) {
    	
	    	curr.bootstrap.results = bootstrap.results[[m]]
	    	
	    	for(i in 1:nrow(curr.bootstrap.results)) {
	    		
	    		curr.adj.matrix = array(0, c(nevents(reconstruction),nevents(reconstruction)))
	    		
	    		for(j in 1:ncol(curr.bootstrap.results)) {
	    			
	    			curr.result = curr.bootstrap.results[i, j]
	    			
	    			for(k in 1:length(curr.result)) {
	    				
	    				curr.val = curr.result[[k]]
	                
	                	for(l in 1:length(curr.val)) {
	                    	if(length(curr.val[l])>1 || curr.val[l] != -1) {
	                        	curr.adj.matrix[curr.val[l], j] = 1
	                    	}
	                	}
	                	
	    			}
	    			
	    		}
	        	
	        	# if I have a perfect match between the reconstructed topologies, increase the count
	        	reconstructed.topology = as.adj.matrix(reconstruction,model=m)[[m]]
	        flag = TRUE;
	        for (j in 1:nrow(reconstructed.topology)) {
		        	for (k in 1:ncol(reconstructed.topology)) {
		        		if(reconstructed.topology[j,k]!=curr.adj.matrix[j,k]) {
		        			flag = FALSE;
		        			next();
		        		}
		        	}
	        }
	        	if(flag==TRUE) {
	            	overall.confidence[[m]] = overall.confidence[[m]] + 1
	            	overall.frequency[[m]] = overall.confidence[[m]] / nboot
	        	}
	    		
	    	}
    	
    }
    
    # save the edge confidence and the frequency of the bootstrap adj.matrix
    for (m in names(as.models(reconstruction))) {
    		
    		curr.adj.matrix = as.adj.matrix(reconstruction,model=m)[[m]];
    	
	    	# save the edge confidence
	    	curr.bootstrap.matrix = bootstrap.adj.matrix[[m]][-1,-1];
	    	curr.edge.confidence = array(0,c(ncol(curr.bootstrap.matrix),nrow(curr.bootstrap.matrix)))
	    	colnames(curr.edge.confidence) = colnames(curr.bootstrap.matrix);
	    	rownames(curr.edge.confidence) = rownames(curr.bootstrap.matrix);
	    	for (i in 1:ncol(curr.bootstrap.matrix)) {
	    			for (j in 1:nrow(curr.bootstrap.matrix)) {
	    				curr.edge.confidence[i,j] = (curr.adj.matrix[i,j]*as.numeric(curr.bootstrap.matrix[i,j]))/nboot
	    			}
	    		}
	    	bootstrap.edge.confidence[[m]] = curr.edge.confidence
	    		
	    		# save the frequency of the bootstrap adj.matrix
	    	curr.bootstrap.matrix = bootstrap.adj.matrix[[m]];
	    	curr.adj.matrix.frequency = array(0,c(ncol(curr.bootstrap.matrix),nrow(curr.bootstrap.matrix)))
	    	colnames(curr.adj.matrix.frequency) = colnames(curr.bootstrap.matrix);
	    	rownames(curr.adj.matrix.frequency) = rownames(curr.bootstrap.matrix);
	   		for (i in 1:ncol(curr.bootstrap.matrix)) {
	   			for (j in 1:nrow(curr.bootstrap.matrix)) {
	    			curr.adj.matrix.frequency[i,j] = as.numeric(as.numeric(curr.bootstrap.matrix[i,j]))/nboot
	    		}
	    	}
    		bootstrap.adj.matrix.frequency[[m]] = curr.adj.matrix.frequency
    
    }
    
    # save the statistics of the bootstrap
    bootstrap.statistics = list(bootstrap.results = bootstrap.results,
    							bootstrap.adj.matrix = list(count = bootstrap.adj.matrix, frequency = bootstrap.adj.matrix.frequency),
    							bootstrap.edge.confidence = bootstrap.edge.confidence,
    							overall.confidence = list(count = overall.confidence, frequency = overall.frequency),
                            bootstrap.settings = list(type = command, nboot = nboot))

    return(bootstrap.statistics)
    
}

#### end of file -- bootstrap.caprese.R


#### decimal.to.binary.tree.R
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
"decimal.to.binary.tree" <-
function(num.decimal, num.bits) {
    #structure where to save the result
    num.binary = rep(0,num.bits);
    #convert the integer decimal number to binary
    pos = 0;
    while(num.decimal>0) {
        #compute the value of the current step
        num.binary[num.bits-pos] = num.decimal %% 2;
        #divide the number by 2 for the next iteration
        num.decimal = num.decimal %/% 2;
        pos = pos + 1;
    }
    return(num.binary);
}

#### end of file -- decimal.to.binary.tree.R


#### estimate.tree.samples.R
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
#error.rates: error rates for false positives and false negatives
#RETURN:
#probabilities: probability of each sample
"estimate.tree.samples" <- function(dataset, reconstructed.topology, estimated.marginal.probabilities, estimated.conditional.probabilities, error.rates) {
    #structure where to save the probabilities of the samples
    probabilities = array(-1,c(nrow(dataset),1));
    #topological properties:
    #1. tree number
    #2. parent
    #3. level in the tree
    topology.structure = array(0,c(nrow(reconstructed.topology),3));
    #go through the subtrees within the topology of four 
    tree.count = 0;
    for (i in 1:nrow(reconstructed.topology)) {
        #if node i has no parents, it is a root
        if(length(which(reconstructed.topology[,i]==1))==0) {
            tree.count = tree.count + 1;
            level = 1;
            #set the parameters for the root
            topology.structure[i,1] = tree.count;
            topology.structure[i,2] = -1;
            topology.structure[i,3] = level;
            curr.node = i;
            #go through this tree
            while (length(curr.node)>0) {
                #move to the next level
                level = level + 1;
                new.node = vector();
                for (j in 1:length(curr.node)) {
                    curr.new.node = which(reconstructed.topology[curr.node[j],]==1);
                    if(length(curr.new.node)>0) {
                        new.node = c(new.node,curr.new.node);
                        for (k in 1:length(curr.new.node)) {
                            #number of the current subtree
                            topology.structure[curr.new.node[k],1] = tree.count;
                            #parent of the current node
                            topology.structure[curr.new.node[k],2] = curr.node[j];
                            #level of this node
                            topology.structure[curr.new.node[k],3] = level;
                        }
                    }
                }
                curr.node = new.node;
            }
        }
    }
    #go through the dataset and evalutate the probability of each sample
    for (i in 1:nrow(dataset)) {
        sample.probability = 1;
        for (j in 1:tree.count) {
            #probability of this subtree (without any knowledge, I set it to 1)
            curr.sample.probability = 1;
            #entries referring to this subtree
            curr.entry = which(topology.structure[,1]==j);
            #samples of each element of this subtree
            curr.sample = dataset[i,curr.entry];
            #parents of each element of this subtree
            curr.parents = topology.structure[curr.entry,2];
            #level of each element of this subtree
            curr.levels = topology.structure[curr.entry,3];
            #set the probability as the one of the root of this tree
            curr.sample.probability = curr.sample.probability * estimated.marginal.probabilities[curr.entry[which(curr.levels==1,arr.ind=TRUE)],1];
            #set the maximum level of this subtree
            max.level = curr.levels[which.max(curr.levels)];
            #if I have at least one event in this sample
            if(length(curr.sample[curr.sample==1])>0) {
                #visit the nodes starting from the lower level
                is.valid = TRUE;
                for (k in max.level:1) {
                    curr.level.nodes = which(curr.levels==k,arr.ind=TRUE);
                    #if I'm not on a root
                    if(k>1) {
                        curr.level.samples = curr.sample[curr.level.nodes];
                        #if I have at least one event at this level
                        if(length(curr.level.samples[curr.level.samples==1])>0) {
                            #I can not have a child without its parent
                            curr.level.parent = curr.parents[curr.level.nodes];
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i,curr.level.parent[p]]==0 && dataset[i,curr.entry[curr.level.nodes[p]]]==1) {
                                    is.valid = FALSE;
                                    break;
                                }
                            }
                        }
                        #if the sample is valid
                        if(is.valid==TRUE) {
                            #add the probability of each edge
                            curr.level.parent = curr.parents[curr.level.nodes];
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i,curr.level.parent[p]]==1 && dataset[i,curr.entry[curr.level.nodes[p]]]==0) {
                                    curr.sample.probability = curr.sample.probability * (1 - estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]],1]);
                                }
                                else if(dataset[i,curr.level.parent[p]]==1 && dataset[i,curr.entry[curr.level.nodes[p]]]==1) {
                                    curr.sample.probability = curr.sample.probability * estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]],1];
                                }
                            }
                        }
                    }
                    if(is.valid==FALSE) {
                        curr.sample.probability = 0;
                        break;
                    }
                }
                if(is.valid==FALSE) {
                    sample.probability = 0;
                    break;
                }
            }
            #if this sample has no events for this tree
            else {
                curr.sample.probability = 1 - curr.sample.probability;
            }
            #update the probability of the topology with the one of this sample
            sample.probability = sample.probability * curr.sample.probability;
            if(sample.probability==0) {
                break;
            }
        }
        probabilities[i,1] = sample.probability;
    }
    #correct the estimation by the error rates
    errors.matrix <- array(0,c(nrow(probabilities),nrow(dataset)));
    for (i in 1:nrow(probabilities)) {
        for (j in 1:nrow(dataset)) {
            curr.sample.x = as.numeric(dataset[i,]);
            curr.sample.y = as.numeric(dataset[j,]);
            errors.matrix[i,j] = (1-error.rates$error.fp)^((1-curr.sample.x)%*%(1-curr.sample.y))*error.rates$error.fp^((1-curr.sample.x)%*%curr.sample.y)*(1-error.rates$error.fn)^(curr.sample.x%*%curr.sample.y)*error.rates$error.fn^(curr.sample.x%*%(1-curr.sample.y));
        }
    }
    probabilities[,1] = as.numeric(as.vector(probabilities)%*%errors.matrix);
    return(probabilities);
}

#### end of file -- estimate.tree.samples.R
