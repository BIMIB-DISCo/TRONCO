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
"bootstrap.caprese" <-
function(dataset, lambda, reconstructed.topology, command = "non-parametric", estimated.marginal.probabilities, estimated.conditional.probabilities, error.rates, nboot) {
    #structure to save the statistics of the bootstrap
    bootstrap.adj.matrix = array(0,c(ncol(dataset)+1,ncol(dataset)+1));
    colnames(bootstrap.adj.matrix) = c("None",colnames(dataset));
    rownames(bootstrap.adj.matrix) = c("None",colnames(dataset));
    #set the dataset if the bootstrap is parametric
    if(command=="parametric") {
    		#define the possible samples given the current number of events
        possible.strings = 2^ncol(dataset);
        err = "";
  		message = "Too many events! Parametric bootstrastap can not be performed."
        err <- tryCatch(curr.dataset <- suppressWarnings(array(0,c(possible.strings,ncol(dataset)))), error = function(e) err <- message);
        if(toString(err) == message) {
        		stop(err, call. = FALSE);
        }
        for (i in 1:possible.strings) {
            curr.dataset[i,] = decimal.to.binary.tree(i-1,ncol(dataset));
        }
        colnames(curr.dataset) = colnames(dataset);
        #define the samples distribution induced by the topology
        samples.probabilities = estimate.tree.samples(curr.dataset,reconstructed.topology,estimated.marginal.probabilities,estimated.conditional.probabilities,error.rates);
    }
    #structure to save the results of the bootstrap
    bootstrap.results = array(-1,c(nboot,ncol(dataset)));
    colnames(bootstrap.results) = colnames(dataset);    
	#perform nboot bootstrap resampling
  
	# create a progress bar
	flush.console()

	pb <- txtProgressBar(1, nboot, style = 3);
  
    for (num in 1:nboot) {
      setTxtProgressBar(pb, num)
    		#performed the bootstrapping procedure
    		if(command=="non-parametric") {
    			#perform the sampling for the current step of bootstrap
    			samples <- sample(1:nrow(dataset),size=nrow(dataset),replace=TRUE);
        		#perform the reconstruction on the bootstrapped dataset
        		check.data = check.dataset(dataset[samples,],FALSE);
    		}
    		else if(command=="parametric") {
    			#perform the sampling for the current step of bootstrap
    			samples <- suppressWarnings(sample(1:nrow(curr.dataset),size=nrow(dataset),replace=TRUE,prob=samples.probabilities));
        		#perform the reconstruction on the bootstrapped dataset
        		check.data = check.dataset(curr.dataset[samples,],FALSE);
    		}
        #if the reconstruction was performed without errors
        if(check.data$is.valid==TRUE) {
        		bootstrapped.dataset = check.data$dataset;
            bootstrapped.topology = caprese.fit(bootstrapped.dataset,lambda,FALSE);
            #set the reconstructed causal edges
            parents.pos = array(-1,c(ncol(bootstrapped.topology$data),1));
            for(i in 1:ncol(bootstrapped.topology$data)) {
            		for(j in 1:ncol(bootstrapped.topology$data)) {
                		if(i!=j && bootstrapped.topology$adj.matrix[i,j]==1) {
						parents.pos[j,1] = i;
                    }
                }
            }
            #get the matched edge in the reconstruction
            matched.idx = match(colnames(bootstrapped.topology$data),colnames(bootstrap.results));
            #if an event has no match, it means it has been merged and I discard it
            parents.pos = parents.pos[!is.na(matched.idx)];
            matched.idx = matched.idx[!is.na(matched.idx)];
            #save the results
            bootstrap.results[num,matched.idx] = parents.pos;
        }
    }
  
	# close progress bar
	close(pb);
  
    #set the statistics of the bootstrap
    for(i in 1:ncol(bootstrap.adj.matrix)) {
        for(j in 1:ncol(bootstrap.adj.matrix)) {
            #if the edge is valid (no self cause)
            if(i!=j) {
                if(i==1 || j==1) {
                    if(j>1) {
                        curr.result = table(bootstrap.results[,j-1]);
                        curr.result = curr.result[names(curr.result)==-1];
                        if(length(curr.result)>0) {
                            bootstrap.adj.matrix[i,j] = curr.result;
                        }
                    }
                }
                else {
                    curr.result = table(bootstrap.results[,j-1]);
                    curr.result = curr.result[names(curr.result)==(i-1)];
                    if(length(curr.result)>0) {
                        bootstrap.adj.matrix[i,j] = curr.result;
                    }
                }
            }
        }
    }
    #set the parent list of the topology previously reconstructed without any bootstrap
    reconstructed.parents = '';
    for(i in 1:ncol(reconstructed.topology)) {
        curr.parent = which(reconstructed.topology[,i]==1);
        if(length(curr.parent)==0) {
            curr.parent = -1;
        }
        reconstructed.parents = paste(reconstructed.parents,toString(curr.parent),sep='');
    }
    #evalutate the overall confidence
    overall.confidence = 0;
    for(i in 1:nrow(bootstrap.results)) {
        reconstructed.boot = '';
        for(j in 1:ncol(bootstrap.results)) {
            reconstructed.boot = paste(reconstructed.boot,toString(bootstrap.results[i,j]),sep='');
        }
        if(reconstructed.parents==reconstructed.boot) {
            overall.confidence = overall.confidence + 1;
        }
    }
    #save the edge confidence
    edge.confidence = (reconstructed.topology*bootstrap.adj.matrix[-1,-1])/nboot;
    #save the confidence from the bootstrap
    confidence = list(overall.value=overall.confidence,overall.frequency=overall.confidence/nboot,bootstrap.values=bootstrap.adj.matrix[,-1],bootstrap.frequencies=bootstrap.adj.matrix[,-1]/nboot);
    #save the settings of the bootstrap
    bootstrap.settings = list(type=command,nboot=nboot);
    #structure to save the results
    bootstrap.statistics = list(reconstructed.topology=reconstructed.topology,confidence=confidence,edge.confidence=edge.confidence,bootstrap.settings=bootstrap.settings);
    return(bootstrap.statistics);
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
