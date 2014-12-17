#### lifting.operators.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Utility function to add the hypotheses
"aux.log" <-
function( ... ) {
	if(exists("settings") && length(settings$data.values$dataset)>0 && length(list(...))>0) {
		clauses = list(...);
		ret = array(0,c(nrow(settings$data.values$dataset),length(clauses)));
		for (i in 1:length(clauses)) {
			#if the clause is given by name, get the column from settings$data.values
			if(typeof(clauses[[i]])=="character") {
				#if I have the label, get the column in the dataset for this event
				col.num = emap(clauses[[i]]);
				if(col.num==-1) {
					stop(paste("Event ",clauses[[i]]," does not exist! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
				else {
					curr.col = settings$data.values$dataset[,col.num];
					if(length(settings$hypotheses$llist)==0) {
						settings$hypotheses$llist = clauses[i];
					}
					else {
						settings$hypotheses$llist = c(settings$hypotheses$llist,clauses[i]);
					}
					if(length(settings$data.values$cardinality)==0) {
						settings$hypotheses$cardinality = cbind(rep.int(1,settings$data.values$num.events),1:settings$data.values$num.events);
						settings$hypotheses$cardinality = rbind(settings$hypotheses$cardinality,c(length(clauses),nrow(settings$hypotheses$cardinality)+1));
					}
					else {
						settings$hypotheses$cardinality = rbind(settings$data.values$cardinality,c(length(clauses),nrow(settings$data.values$cardinality)+1));
					}
					assign("settings",settings,envir=.GlobalEnv);
				}
			}
			#otherwise I already have the column as a vector
			else {
				#if it is a vector
				curr.col = clauses[[i]];
			}
			ret[,i] = curr.col;
		}
		return(ret);
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
	return(NA);
}

# XOR hypothesis
"XOR" <-
function( ... ) {
	if(exists("settings") && length(settings$data.values$dataset)>0 && length(list(...))>0) {
		#get the vector of the clauses of the formula from settings$data.values
		ret = aux.log(...);
		#evaluate the XOR operator
		result = rep(0,nrow(settings$data.values$dataset));
		for (i in 1:nrow(settings$data.values$dataset)) {
			result[i] = sum(ret[i,]);
			if(result[i]>1) {
				result[i] = 0;
			}
		}
		return(as.integer(result));
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
	return(NA);
}

# OR hypothesis
"OR" <-
function( ... ) {
	if(exists("settings") && length(settings$data.values$dataset)>0 && length(list(...))>0) {
		#get the vector of the clauses of the formula from settings$data.values
		ret = aux.log(...);
		#evaluate the OR operator
		result = rep(0,nrow(settings$data.values$dataset));
		for (i in 1:nrow(settings$data.values$dataset)) {
			result[i] = sum(ret[i,]);
			if(result[i]>1) {
				result[i] = 1;
			}
		}
		return(as.integer(result));
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
	return(NA);
}

# AND hypothesis
"AND" <-
function( ... ) {
	if(exists("settings") && length(settings$data.values$dataset)>0 && length(list(...))>0) {
		#get the vector of the clauses of the formula from settings$data.values
		ret = aux.log(...);
		#evaluate the AND operator
		result = rep(0,nrow(settings$data.values$dataset));
		for (i in 1:nrow(settings$data.values$dataset)) {
			result[i] = sum(ret[i,]);
			if(result[i]<ncol(ret)) {
				result[i] = 0;
			}
			else {
				result[i] = 1;
			}
		}
		return(as.integer(result));
	}
	else {
		stop("Either the dataset or the formula is not provided! No hypothesis will be created.");
	}
	return(NA);
}

#### end of file -- lifting.operators.R
