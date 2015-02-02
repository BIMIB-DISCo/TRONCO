#### bootstrap.capri.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#perform non-parametric or parametric bootstrap to evalutate the confidence of the reconstruction
#INPUT:
#dataset: a dataset describing a progressive phenomenon
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
"bootstrap.capri" <-
function(dataset, command.capri, do.boot, nboot.capri, pvalue, reconstructed.topology.pf, reconstructed.topology.bic, command = "non-parametric", estimated.marginal.probabilities.pf, estimated.conditional.probabilities.pf, parents.pos.pf, error.rates.pf, estimated.marginal.probabilities.bic, estimated.conditional.probabilities.bic, parents.pos.bic, error.rates.bic, nboot) {
    #structure to save the statistics of the bootstrap
    bootstrap.adj.matrix.pf = array(0,c(ncol(dataset)+1,ncol(dataset)+1));
    colnames(bootstrap.adj.matrix.pf) = c("None",colnames(dataset));
    rownames(bootstrap.adj.matrix.pf) = c("None",colnames(dataset));
    bootstrap.adj.matrix.bic = array(0,c(ncol(dataset)+1,ncol(dataset)+1));
    colnames(bootstrap.adj.matrix.bic) = c("None",colnames(dataset));
    rownames(bootstrap.adj.matrix.bic) = c("None",colnames(dataset));
    #structure to save the results of the bootstrap
    bootstrap.results.pf = array(list(-1),c(nboot,ncol(dataset)));
    colnames(bootstrap.results.pf) = colnames(dataset);
    bootstrap.results.bic = array(list(-1),c(nboot,ncol(dataset)));
    colnames(bootstrap.results.bic) = colnames(dataset);
	#perform nboot bootstrap resampling
  
  	# create a progress bar
	pb <- txtProgressBar(1, nboot, style = 3);
  
    for (num in 1:nboot) {
      setTxtProgressBar(pb, num)
		#performed the bootstrapping procedure
		if(command=="non-parametric") {
			#perform the sampling for the current step of bootstrap
			samples <- sample(1:nrow(dataset),size=nrow(dataset),replace=TRUE);
			#perform the reconstruction on the bootstrapped dataset
			check.data = check.dataset(dataset[samples,],FALSE);
			#if the reconstruction was performed without errors
			if(check.data$is.valid==TRUE) {
				bootstrapped.dataset = check.data$dataset;
				bootstrapped.hypotheses = check.data$hypotheses;
				bootstrapped.topology = capri.fit(bootstrapped.dataset,bootstrapped.hypotheses,command.capri,do.boot,nboot.capri,pvalue,FALSE);
				#set the reconstructed causal edges
				parents.pos.pf = array(list(),c(ncol(bootstrapped.topology$data),1));
				parents.pos.bic = array(list(),c(ncol(bootstrapped.topology$data),1));
				for(i in 1:ncol(bootstrapped.topology$data)) {
					for(j in 1:ncol(bootstrapped.topology$data)) {
						if(i!=j && bootstrapped.topology$adj.matrix$adj.matrix.pf[i,j]==1) {
							parents.pos.pf[j,1] = list(c(unlist(parents.pos.pf[j,1]),i));
						}
						if(i!=j && bootstrapped.topology$adj.matrix$adj.matrix.bic[i,j]==1) {
							parents.pos.bic[j,1] = list(c(unlist(parents.pos.bic[j,1]),i));
						}
					}
				}
				parents.pos.pf[unlist(lapply(parents.pos.pf,is.null))] = list(-1);
				parents.pos.bic[unlist(lapply(parents.pos.bic,is.null))] = list(-1);
				#get the matched edge in the reconstruction
				matched.idx.pf = match(colnames(bootstrapped.topology$data),colnames(bootstrap.results.pf));
				matched.idx.bic = match(colnames(bootstrapped.topology$data),colnames(bootstrap.results.bic));
				#if an event has no match, it means it has been merged and I discard it
				parents.pos.pf = parents.pos.pf[!is.na(matched.idx.pf)];
				parents.pos.bic = parents.pos.bic[!is.na(matched.idx.bic)];
				matched.idx.pf = matched.idx.pf[!is.na(matched.idx.pf)];
				matched.idx.bic = matched.idx.bic[!is.na(matched.idx.bic)];
				#save the results
				bootstrap.results.pf[num,matched.idx.pf] = parents.pos.pf;
				bootstrap.results.bic[num,matched.idx.bic] = parents.pos.bic;
			}
		}
		else if(command=="parametric") {
			if(num==1) {
				#define the possible samples given the current number of events
				possible.strings = 2^ncol(dataset);
				err = "";
				message = "Too many events! Parametric bootstrastap can not be performed."
				err <- tryCatch(curr.dataset <- suppressWarnings(array(0,c(possible.strings,ncol(dataset)))), error = function(e) err <- message);
				if(toString(err) == message) {
					stop(err, call. = FALSE);
				}
				for (i in 1:possible.strings) {
					curr.dataset[i,] = decimal.to.binary.dag(i-1,ncol(dataset));
				}
				colnames(curr.dataset) = colnames(dataset);
				#define the samples distribution induced by the topology
				samples.probabilities.pf = estimate.dag.samples(curr.dataset,reconstructed.topology.pf,estimated.marginal.probabilities.pf,estimated.conditional.probabilities.pf,parents.pos.pf,error.rates.pf);
				samples.probabilities.bic = estimate.dag.samples(curr.dataset,reconstructed.topology.bic,estimated.marginal.probabilities.bic,estimated.conditional.probabilities.bic,parents.pos.bic,error.rates.bic);
			}
			#perform the sampling for the current step of bootstrap
			samples.pf <- suppressWarnings(sample(1:nrow(curr.dataset),size=nrow(dataset),replace=TRUE,prob=samples.probabilities.pf));
			samples.bic <- suppressWarnings(sample(1:nrow(curr.dataset),size=nrow(dataset),replace=TRUE,prob=samples.probabilities.bic));
			#perform the reconstruction on the bootstrapped dataset
			check.data.pf = check.dataset(curr.dataset[samples.pf,],FALSE);
			check.data.bic = check.dataset(curr.dataset[samples.bic,],FALSE);
			#if the reconstruction was performed without errors for the prima facie topology
			if(check.data.pf$is.valid==TRUE) {
				bootstrapped.dataset = check.data.pf$dataset;
				bootstrapped.hypotheses = check.data.pf$hypotheses;
				bootstrapped.topology = capri.fit(bootstrapped.dataset,bootstrapped.hypotheses,command.capri,do.boot,nboot.capri,pvalue,FALSE);
				#set the reconstructed causal edges
				parents.pos.pf = array(list(),c(ncol(bootstrapped.topology$data),1));
				for(i in 1:ncol(bootstrapped.topology$data)) {
					for(j in 1:ncol(bootstrapped.topology$data)) {
						if(i!=j && bootstrapped.topology$adj.matrix$adj.matrix.pf[i,j]==1) {
							parents.pos.pf[j,1] = list(c(unlist(parents.pos.pf[j,1]),i));
						}
					}
				}
				parents.pos.pf[unlist(lapply(parents.pos.pf,is.null))] = list(-1);
				#get the matched edge in the reconstruction
				matched.idx.pf = match(colnames(bootstrapped.topology$data),colnames(bootstrap.results.pf));
				#if an event has no match, it means it has been merged and I discard it
				parents.pos.pf = parents.pos.pf[!is.na(matched.idx.pf)];
				matched.idx.pf = matched.idx.pf[!is.na(matched.idx.pf)];
				#save the results
				bootstrap.results.pf[num,matched.idx.pf] = parents.pos.pf;
			}
			#if the reconstruction was performed without errors for the causal topology
			if(check.data.bic$is.valid==TRUE) {
				bootstrapped.dataset = check.data.bic$dataset;
				bootstrapped.hypotheses = check.data.bic$hypotheses;
				bootstrapped.topology = capri.fit(bootstrapped.dataset,bootstrapped.hypotheses,command.capri,do.boot,nboot.capri,pvalue,FALSE);
				#set the reconstructed causal edges
				parents.pos.bic = array(list(),c(ncol(bootstrapped.topology$data),1));
				for(i in 1:ncol(bootstrapped.topology$data)) {
					for(j in 1:ncol(bootstrapped.topology$data)) {
						if(i!=j && bootstrapped.topology$adj.matrix$adj.matrix.bic[i,j]==1) {
							parents.pos.bic[j,1] = list(c(unlist(parents.pos.bic[j,1]),i));
						}
					}
				}
				parents.pos.bic[unlist(lapply(parents.pos.bic,is.null))] = list(-1);
				#get the matched edge in the reconstruction
				matched.idx.bic = match(colnames(bootstrapped.topology$data),colnames(bootstrap.results.bic));
				#if an event has no match, it means it has been merged and I discard it
				parents.pos.bic = parents.pos.bic[!is.na(matched.idx.bic)];
				matched.idx.bic = matched.idx.bic[!is.na(matched.idx.bic)];
				#save the results
				bootstrap.results.bic[num,matched.idx.bic] = parents.pos.bic;
			}
		}
    }
  
    # close progress bar
    close(pb);
  
    #set the statistics of the bootstrap
    for(i in 2:ncol(bootstrap.adj.matrix.pf)) {
		curr.result.pf = bootstrap.results.pf[,i-1];
        curr.result.bic = bootstrap.results.bic[,i-1];
        for(j in 1:length(curr.result.pf)) {
			curr.pf = curr.result.pf[[j]];
            for(k in 1:length(curr.pf)) {
				if(curr.pf[k]==-1) {
					bootstrap.adj.matrix.pf[1,i] = bootstrap.adj.matrix.pf[1,i] + 1;
				}
				else {
					bootstrap.adj.matrix.pf[curr.pf[k]+1,i] = bootstrap.adj.matrix.pf[curr.pf[k]+1,i] + 1;
				}
            }
            curr.bic = curr.result.bic[[j]];
            for(k in 1:length(curr.bic)) {
				if(curr.bic[k]==-1) {
					bootstrap.adj.matrix.bic[1,i] = bootstrap.adj.matrix.bic[1,i] + 1;
				}
				else {
					bootstrap.adj.matrix.bic[curr.bic[k]+1,i] = bootstrap.adj.matrix.bic[curr.bic[k]+1,i] + 1;
				}
            }
        }
    }
    #evalutate the overall confidence
    overall.confidence.pf = 0;
    overall.confidence.bic = 0;
    cont = 0;
    for(i in 1:nrow(bootstrap.results.pf)) {
		curr.adj.matrix.pf = array(0,c(ncol(dataset),ncol(dataset)));
		curr.adj.matrix.bic = array(0,c(ncol(dataset),ncol(dataset)));
		for(j in 1:ncol(bootstrap.results.pf)) {
			curr.result.pf = bootstrap.results.pf[i,j];
			curr.result.bic = bootstrap.results.bic[i,j];
			for(k in 1:length(curr.result.pf)) {
				curr.pf = curr.result.pf[[k]];
				for(l in 1:length(curr.pf)) {
					if(curr.pf[l]!=-1) {
						curr.adj.matrix.pf[curr.pf[l],j] = 1;
					}
				}
				curr.bic = curr.result.bic[[k]];
				for(l in 1:length(curr.bic)) {
					if(curr.bic[l]!=-1) {
						curr.adj.matrix.bic[curr.bic[l],j] = 1;
					}
				}
			}
		}
		if(sum(reconstructed.topology.pf-curr.adj.matrix.pf)==0) {
			overall.confidence.pf = overall.confidence.pf + 1;
		}
		if(sum(reconstructed.topology.bic-curr.adj.matrix.bic)==0) {
			overall.confidence.bic = overall.confidence.bic + 1;
		}
    }
    #save the reconstructed topologies
    reconstructed.topology = list(reconstructed.topology.pf=reconstructed.topology.pf,reconstructed.topology.bic=reconstructed.topology.bic);
    #save the edge confidence
    edge.confidence.pf = (reconstructed.topology.pf*bootstrap.adj.matrix.pf[-1,-1])/nboot;
    edge.confidence.bic = (reconstructed.topology.bic*bootstrap.adj.matrix.bic[-1,-1])/nboot;
    edge.confidence = list(edge.confidence.pf=edge.confidence.pf,edge.confidence.bic=edge.confidence.bic);
    #save the confidence from the bootstrap
    confidence.pf = list(overall.value.pf=overall.confidence.pf,overall.frequency.pf=overall.confidence.pf/nboot,bootstrap.values.pf=bootstrap.adj.matrix.pf[,-1],bootstrap.frequencies.pf=bootstrap.adj.matrix.pf[,-1]/nboot);
    confidence.bic = list(overall.value.bic=overall.confidence.bic,overall.frequency.bic=overall.confidence.bic/nboot,bootstrap.values.bic=bootstrap.adj.matrix.bic[,-1],bootstrap.frequencies.bic=bootstrap.adj.matrix.bic[,-1]/nboot);
    confidence = list(confidence.pf=confidence.pf,confidence.bic=confidence.bic);    
    #save the settings of the bootstrap
    bootstrap.settings = list(type=command,nboot=nboot);
    #structure to save the results
    bootstrap.statistics = list(reconstructed.topology=reconstructed.topology,confidence=confidence,edge.confidence=edge.confidence,bootstrap.settings=bootstrap.settings);
    return(bootstrap.statistics);
}

#### end of file -- bootstrap.capri.R
