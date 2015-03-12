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
