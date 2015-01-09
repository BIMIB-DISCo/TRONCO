#### hypothesis.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Add a new hypothesis by creating a new causal event and adding it to the dateset
"hypothesis.add" <-
function( label.formula, lifted.formula, label.effect  ) {
	if(exists("settings") && length(settings$data.values)>0 && length(settings$visualization)>0 && length(settings$hypotheses)>0) {
		#this is required because R performs lazy call by name
		lifted.formula;
		#look for duplicated causes in the formula
		if(anyDuplicated(settings$hypotheses$llist)>0) {
			settings$hypotheses$llist = list();
			settings$hypotheses$cardinality = NULL;
			assign("settings",settings,envir=.GlobalEnv);
			stop(paste("There are duplicated causes in the formula! No hypothesis will be created.",sep=''));
		}
		#* is a special label.effect which indicates to use all the events as effects for the formula
		if(length(label.effect)==1 && label.effect=="*") {
			label.effect = names(settings$data.values$dataset)[1:(length(names(settings$data.values$dataset))-settings$hypotheses$num.hypotheses)];
			if(length(settings$hypotheses$llist)==0) {
				settings$hypotheses$cardinality = NULL;
				assign("settings",settings,envir=.GlobalEnv);
				stop(paste("No valid clause has been provided in the formula! No hypothesis will be created.",sep=''));
			}
			else {
				#the events can not be both causes and effects for the formula to be well-formed
				for(i in 1:length(settings$hypotheses$llist)) {
					val = which(label.effect==settings$hypotheses$llist[[i]]);
					label.effect = label.effect[-val];
				}
			}
			if(length(label.effect)==0) {
				settings$hypotheses$llist = list();
				settings$hypotheses$cardinality = NULL;
				assign("settings",settings,envir=.GlobalEnv);
				stop(paste("No valid effect provided in the formula! No hypothesis will be created.",sep=''));
			}
		}
		#check the formula to be well-formed
		all.col.nums = vector();
		if(length(label.effect)==0) {
			settings$hypotheses$llist = list();
			settings$hypotheses$cardinality = NULL;
			assign("settings",settings,envir=.GlobalEnv);
			stop(paste("No valid effect provided in the formula! No hypothesis will be created.",sep=''));
		}
		else {
			#check the effects of the formula to be well-formed
			for (i in 1:length(label.effect)) {
				col.num = emap(label.effect[[i]]);
				all.col.nums[length(all.col.nums)+1] = col.num;
				#check the effect to be a valid event
				if(col.num==-1) {
					settings$hypotheses$llist = list();
					settings$hypotheses$cardinality = NULL;
					assign("settings",settings,envir=.GlobalEnv);
					stop(paste("Event ",label.effect[[i]]," does not exist! The effect is bad formed and no hypothesis will be created.",sep=''));
				}
				#check the formula to be well-formed
				#if the effect is in the formula, the formula is not well-formed
				for (j in 1:length(settings$hypotheses$llist)) {
					if(label.effect[[i]]==settings$hypotheses$llist[[j]]) {
						settings$hypotheses$llist = list();
						settings$hypotheses$cardinality = NULL;
						assign("settings",settings,envir=.GlobalEnv);
						stop(paste("The effect is in the formula! The formula is bad formed and no hypothesis will be created.",sep=''));
					}
				}
			}
		}
		#look for duplicated effects in the formula
		if(anyDuplicated(all.col.nums)>0) {
			settings$hypotheses$llist = list();
			settings$hypotheses$cardinality = NULL;
			assign("settings",settings,envir=.GlobalEnv);
			stop(paste("There are duplicated effects in the formula! No hypothesis will be created.",sep=''));
		}
		#check that the we are not duplicating any name by adding the new hypothesis
		if(length(which(names(settings$data.values$dataset)==label.formula))>0) {
			settings$hypotheses$llist = list();
			settings$hypotheses$cardinality = NULL;
			assign("settings",settings,envir=.GlobalEnv);
			stop(paste("Hypothesis ",label.formula," already exists! The hypothesis will not be created.",sep=''));
		}
		#add the hypothesis to the dataset
		tmp.data.values = settings$data.values$dataset;
		settings$data.values$dataset = cbind(settings$data.values$dataset,lifted.formula);
		assign("settings",settings,envir=.GlobalEnv);
		#check that the formula is valid
		#structure to compute the observed and observed joint probabilities
		pair.count <- array(0, dim=c(ncol(settings$data.values$dataset),ncol(settings$data.values$dataset)));
		#compute the probabilities on the dataset
		for(i in 1:ncol(settings$data.values$dataset)) {
			for(j in 1:ncol(settings$data.values$dataset)) {
				val1 = settings$data.values$dataset[ ,i];
				val2 = settings$data.values$dataset[ ,j];
				pair.count[i,j] = (t(val1) %*% val2);
			}
		}
		#marginal.probs is an array of the observed marginal probabilities
		marginal.probs <- array(as.matrix(diag(pair.count)/nrow(settings$data.values$dataset)),dim=c(ncol(settings$data.values$dataset),1));
		#joint.probs is an array of the observed joint probabilities
		joint.probs <- as.matrix(pair.count/nrow(settings$data.values$dataset));
		#check that the probability of the formula is in (0,1)
		if(marginal.probs[ncol(settings$data.values$dataset)]==0 || marginal.probs[ncol(settings$data.values$dataset)]==1) {
			settings$data.values$dataset = tmp.data.values;
			assign("settings",settings,envir=.GlobalEnv);
			stop(paste("The marginal probability of the formula is not strictly in (0,1)! The formula is bad formed and no hypothesis will be created.",sep=''));
		}
		#check that the formula does not duplicate any existing column
		i = ncol(settings$data.values$dataset);
		for(j in 1:ncol(settings$data.values$dataset)) {
			#if the edge is valid, i.e., not self cause
			if(i!=j) {
				#if the two considered events are not distinguishable
				if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) {
					settings$data.values$dataset = tmp.data.values;
					assign("settings",settings,envir=.GlobalEnv);
					stop(paste("The formula duplicates an already existing event! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
			}
		}
		#now I can finally add the hypothesis
		names(settings$data.values$dataset)[length(settings$data.values$dataset)] = label.formula;
		settings$visualization$labels = c(settings$visualization$labels,label.formula);
		tmp = strsplit(label.formula,":")[[1]];
		settings$visualization$visualized.labels = c(settings$visualization$visualized.labels,tmp[1]);
		settings$visualization$colors = c(settings$visualization$colors,"lightyellow");
		settings$hypotheses$num.hypotheses = settings$hypotheses$num.hypotheses + 1;
		settings$hypotheses$llist = list();
		settings$data.values$cardinality = settings$hypotheses$cardinality;
		settings$hypotheses$cardinality = NULL;
		#create the list of added hypotheses
		if(length(settings$hypotheses$hlist)==0) {
			settings$hypotheses$hlist = vector();
		}
		#add the new hypothesis to the list
		for (i in 1:length(label.effect)) {
			col.num = emap(label.effect[[i]]);
			settings$hypotheses$hlist = c(settings$hypotheses$hlist,settings$hypotheses$num.hypotheses,col.num);
		}
		assign("settings",settings,envir=.GlobalEnv);
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
}

#### end of file -- hypothesis.add.R
