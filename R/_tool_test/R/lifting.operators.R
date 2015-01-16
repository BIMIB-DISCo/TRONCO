#### lifting.operators.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Utility function to add the hypotheses
"aux.log" <-
function( dataset, function.name, ... ) {
	if(!is.null(dataset) && length(list(...))>0) {
		clauses = list(...);
		curr_dataset = array(0,c(nrow(dataset),length(clauses)));
		hypotheses = list();
		function.effects = list();
		for (i in 1:length(clauses)) {
			#if the clause is given by name, get the column from the dataset
			if(typeof(clauses[[i]])=="character") {
				#if I have the label, get the column in the dataset for this event
				col.num = emap(clauses[[i]],dataset);
				if(col.num==-1) {
					stop(paste("Event ",clauses[[i]]," does not exist! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
				else {
					curr_dataset[,i] = dataset[,col.num];
					if(length(hypotheses$llist)==0) {
						hypotheses$llist = clauses[i];
					}
					else {
						hypotheses$llist = list(c(unlist(hypotheses$llist),clauses[[i]]));
					}
					function.name = paste(function.name,"_",clauses[[i]],sep="");
					function.effects = c(function.effects,clauses[[i]]);
				}
			}
			#otherwise I already have the column as a vector
			else {
				#if it is a list
				curr_dataset[,i] = clauses[[i]]$formula;
				if(length(hypotheses$llist)==0) {
					hypotheses$llist = clauses[[i]]$hypotheses$llist;
				}
				else {
					hypotheses$llist = list(c(unlist(clauses[[i]]$hypotheses$llist),unlist(hypotheses$llist)));
				}
				function.name = paste(function.name,"_",clauses[[i]]$function.name,sep="");
				function.effects = c(function.effects,clauses[[i]]$function.name);
			}
		}
		result = list(curr_dataset=curr_dataset,hypotheses=hypotheses,function.name=function.name,function.effects=function.effects);
		#save the new edges
		for(k in 1:length(result$function.effects)) {
			lifting.edges = rbind(lifting.edges,c(result$function.name,result$function.effects[[k]]));
			assign("lifting.edges",lifting.edges,envir=.GlobalEnv);
		}
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA);
}

# AND hypothesis
"AND" <-
function( ... ) {
	#look for a global variable named lifting.dataset
	dataset = lifting.dataset;
	if(!is.null(dataset) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,"AND",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.effects = result$function.effects;
		#evaluate the AND operator
		formula = rep(0,nrow(dataset));
		for (i in 1:nrow(dataset)) {
			formula[i] = sum(curr_dataset[i,]);
			if(formula[i]<ncol(curr_dataset)) {
				formula[i] = 0;
			}
			else {
				formula[i] = 1;
			}
		}
		formula = as.integer(formula);
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.effects=function.effects);
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA);
}

# OR hypothesis
"OR" <-
function( ... ) {
	#look for a global variable named lifting.dataset
	dataset = lifting.dataset;
	if(!is.null(dataset) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,"OR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.effects = result$function.effects;
		#evaluate the OR operator
		formula = rep(0,nrow(dataset));
		for (i in 1:nrow(dataset)) {
			formula[i] = sum(curr_dataset[i,]);
			if(formula[i]>1) {
				formula[i] = 1;
			}
		}
		formula = as.integer(formula);
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.effects=function.effects);
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA);
}

# XOR hypothesis
"XOR" <-
function( ... ) {
	#look for a global variable named lifting.dataset
	dataset = lifting.dataset;
	if(!is.null(dataset) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,"XOR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.effects = result$function.effects;
		#evaluate the XOR operator
		formula = rep(0,nrow(dataset));
		for (i in 1:nrow(dataset)) {
			formula[i] = sum(curr_dataset[i,]);
			if(formula[i]>1) {
				formula[i] = 0;
			}
		}
		formula = as.integer(formula);
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.effects=function.effects);
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA);
}

#### end of file -- lifting.operators.R
