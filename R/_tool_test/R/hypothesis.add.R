#### hypothesis.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Add a new hypothesis by creating a new causal event and adding it to the dateset
"hypothesis.add" <-
function( dataset, label.formula, lifted.formula, label.effect, hypotheses = NA ) {
	if(!is.null(dataset)) {
		#get the series of calls of the recursive function to generate the lifted formula
		hstructure = toString(as.list(match.call())$lifted.formula);
		#save the lifted dataset and its hypotheses for the current formula
		curr_formula = lifted.formula$formula;
		curr_hypotheses = lifted.formula$hypotheses;
		if(!is.na(hypotheses[1])) {
			num.hypotheses = hypotheses$num.hypotheses;
		}
		else {
			num.hypotheses = 0;
		}
		#* is a special label.effect which indicates to use all the events as effects for this formula
		if(label.effect[1]=="*") {
			label.effect = names(dataset)[1:(length(names(dataset))-num.hypotheses)];
			#any event can not be both causes and effects for the formula to be well-formed
			label.effect = label.effect[-which((label.effect%in%unlist(curr_hypotheses$llist)))];
			if(length(label.effect)==0) {
				stop(paste("No valid effect provided in the formula! No hypothesis will be created.",sep=''));
			}
		}
		#check the formula to be well-formed
		all.col.nums = vector();
		if(length(label.effect)==0) {
			stop(paste("No valid effect provided in the formula! No hypothesis will be created.",sep=''));
		}
		else {
			#check the effects of the formula to be well-formed
			for (i in 1:length(label.effect)) {
				col.num = emap(label.effect[[i]],dataset);
				#check the effect to be a valid event
				if(col.num==-1) {
					stop(paste("Event ",label.effect[[i]]," does not exist! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
				all.col.nums[length(all.col.nums)+1] = col.num;
				#check the formula to be well-formed
				#if the effect is in the formula, the formula is not well-formed
				if(length(which(unlist(curr_hypotheses$llist)%in%unlist(label.effect[[i]])))>0) {
					stop(paste("The effect is in the formula! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
			}
		}
		#look for duplicated effects in the formula
		if(anyDuplicated(all.col.nums)>0) {
			stop(paste("There are duplicated effects in the formula! No hypothesis will be created.",sep=''));
		}
		#check that the we are not duplicating any name by adding the new hypothesis
		if(length(which(names(dataset)==label.formula))>0) {
			stop(paste("Hypothesis ",label.formula," already exists! No hypothesis will not be created.",sep=''));
		}
		#add the hypothesis to the dataset
		dataset = cbind(dataset,curr_formula);		
		#check that the formula is valid according to Suppes' theory
		#structure to compute the observed and observed joint probabilities
		pair.count <- array(0, dim=c(ncol(dataset),ncol(dataset)));
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
		#check that the probability of the formula is in (0,1)
		if(marginal.probs[ncol(dataset)]==0 || marginal.probs[ncol(dataset)]==1) {
			stop(paste("The marginal probability of the formula is not strictly in (0,1)! The formula is not valid and no hypothesis will be created.",sep=''));
		}
		#check that the formula does not duplicate any existing column
		i = ncol(dataset);
		for(j in 1:ncol(dataset)) {
			#if the edge is valid, i.e., not self cause
			if(i!=j) {
				#if the two considered events are not distinguishable
				if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) {
					stop(paste("The formula duplicates an already existing event! The formula is not valid and no hypothesis will be created.",sep=''));
				}
			}
		}
		#now I can finally add the hypothesis
		names(dataset)[length(dataset)] = label.formula;
		if(is.na(hypotheses[1])) {
			hypotheses = list();
		}
		hypotheses$num.hypotheses = num.hypotheses + 1;
		#create the list of added hypotheses
		if(length(hypotheses$hlist)==0) {
			hypotheses$hlist = vector();
		}
		#add the new hypothesis to the list
		for (i in 1:length(label.effect)) {
			col.num = emap(label.effect[[i]],dataset);
			hypotheses$hlist = rbind(hypotheses$hlist,t(c(ncol(dataset),col.num)));
			if(is.null(colnames(hypotheses$hlist))) {
				colnames(hypotheses$hlist) = c("cause","effect");
			}
		}
		#create the list of hypotheses' structures
		if(length(hypotheses$hstructure)==0) {
			hypotheses$hstructure = list();
		}
		hypotheses$hstructure = c(hypotheses$hstructure,hstructure);
		#return the result as a list
		result = list(dataset=dataset,hypotheses=hypotheses);
		return(result);
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
	return(NA);
}

#### end of file -- hypothesis.add.R
