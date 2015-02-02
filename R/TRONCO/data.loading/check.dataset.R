#### check.dataset.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#check if the dataset is valid accordingly to the probability raising
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#verbose: should I print the warnings? Yes if TRUE, no otherwise
#RETURN:
#valid.dataset: a dataset valid accordingly to the probability raising
"check.dataset" <-
function( dataset, verbose ) {
    #perform the preprocessing only if I have at least two binary events and two samples
    if(length(ncol(dataset))>0 && ncol(dataset)>1 && length(nrow(dataset))>0 && nrow(dataset)>1 && length(dataset[dataset==0|dataset==1])==nrow(dataset)*ncol(dataset)) {
        #structure to compute the observed and observed joint probabilities
        pair.count <- array(0, dim=c(ncol(dataset), ncol(dataset)));
        #compute the probabilities on the dataset
        for(i in 1:ncol(dataset)) {
            for(j in 1:ncol(dataset)) {
                val1 = dataset[ ,i];
                val2 = dataset[ ,j];
                pair.count[i,j] = (t(val1) %*% val2);
            }
        }
        #marginal.probs is an array of the observed marginal probabilities
        marginal.probs <- array(as.matrix(diag(pair.count)/nrow(dataset)),dim=c(ncol(dataset),1));
        #joint.probs is an array of the observed joint probabilities
        joint.probs <- as.matrix(pair.count/nrow(dataset));
        #remove the events with marginal observed probability not strictly in (0,1), if any
        valid.cols = marginal.probs[,1]>0 & marginal.probs[,1]<1;
        not.valid = which(!valid.cols);
        #save the list of the events to be dropped, if any
        removed.events = NA;
        removed.num = 0;
        if(length(not.valid)>0) {
            removed.events = vector();
            for(i in 1:length(not.valid)) {
                removed.num = removed.num + 1;
                removed.events[removed.num] = not.valid[i];
                if(verbose==TRUE) {
                    warning(paste("Event of column ",toString(not.valid[i])," is not valid and will be discarded.",sep=""));
                }
            }
        }
        #merge groups of events if they are not distinguishable
        merged.events = NA;
        merged.num = 0;
        for(i in 1:ncol(dataset)) {
            for(j in i:ncol(dataset)) {
                #if the edge is valid, i.e., not self cause and not in removed.events
                if(i!=j && (removed.num==0 || (any(removed.events==i)==FALSE && any(removed.events==j)==FALSE))) {
                    #if the two considered events are not distinguishable
                    if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) {
                        #add the events to the merge list
                        if(merged.num==0) {
                            merged.events = vector();
                        }
                        merged.num = merged.num + 1;
                        merged.events[merged.num] = i;
                        merged.num = merged.num + 1;
                        merged.events[merged.num] = j;
                        #add the second event to the drop list
                        if(removed.num==0) {
                            removed.events = vector();
                        }
                        removed.num = removed.num + 1;
                        removed.events[removed.num] = j;
                        #merge the names in header of the data frame
                        names(dataset)[i] = paste(names(dataset)[i],"_and_",names(dataset)[j],sep="");
                        names(dataset)[j] = names(dataset)[i];
                        if(verbose==TRUE) {
                            warning(paste("Events of columns ",toString(i)," and ",toString(j)," are not distinguishable and they will be merged.",sep=""));
                        }
                    }
                }
            }
        }
        #save merged.events to an array
        if(merged.num>0) {
            merged.events = t(array(merged.events, dim=c(2,length(merged.events)/2)));
        }
        #save the valid dataset
        if(removed.num>0) {
            valid.dataset = dataset[,-removed.events];
            valid.marginal.probs = marginal.probs[-removed.events,];
            valid.marginal.probs = array(valid.marginal.probs, dim=c(length(valid.marginal.probs),1));
            valid.joint.probs = joint.probs[-removed.events,-removed.events];
        }
        else {
            valid.dataset = dataset;
            valid.marginal.probs = marginal.probs;
            valid.marginal.probs = array(valid.marginal.probs, dim=c(length(valid.marginal.probs),1));
            valid.joint.probs = joint.probs;
        }
        #if at this point I still have at least two events, the dataset is valid
        if(length(ncol(valid.dataset))>0 && ncol(valid.dataset)>1) {
            invalid.events = list(removed.events=removed.events,merged.events=merged.events);
            valid.dataset = list(dataset=valid.dataset,hypotheses=NA,invalid.events=invalid.events,marginal.probs=valid.marginal.probs,joint.probs=valid.joint.probs,is.valid=TRUE);
        }
        #if the dataset is not valid, we stop here
        else {
            valid.dataset = list(dataset=NA,hypotheses=NA,invalid.events=NA,marginal.probs=NA,joint.probs=NA,is.valid=FALSE);
        }
    }
    #if the dataset is not valid, we stop here
    else {
        if(verbose==TRUE) {
            warning("The dataset must contain at least two binary events and two samples.");
        }
        valid.dataset = list(dataset=NA,hypotheses=NA,invalid.events=NA,marginal.probs=NA,joint.probs=NA,is.valid=FALSE);
    }
    return(valid.dataset);
}

#### end of file -- check.dataset.R
