#### lifting.operators.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Utility function to add the hypotheses
"aux.log" <-
function( dataset, annotations, function.name, ... ) {
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		clauses = list(...);
		curr_dataset = array(0,c(nrow(dataset),length(clauses)));
		hypotheses = list();
		function.inputs = list();
		for (i in 1:length(clauses)) {
			#if the clause is given by name, get the column from the dataset
			if(typeof(clauses[[i]])=="character") {
				#if I have the label, get the column in the dataset for this event
				col.num = -1;
				if(length(clauses[[i]])==1) {
					event.map = emap(c(clauses[[i]],"*"),dataset,annotations);
					col.num = event.map$col.num;
				}
				else if(length(clauses[[i]])==2) {
					event.map = emap(clauses[[i]],dataset,annotations);
					col.num = event.map$col.num;
				}
				if(col.num==-1) {
					stop(paste("One or more undefined events have been found! The formula is bad formed and no hypothesis will be created.",sep=''));
				}
				else {
					curr_dataset[,i] = dataset[,col.num[1]];
					if(length(col.num)>1) {
						curr_dataset = cbind(curr_dataset,dataset[,col.num[2:length(col.num)]]);
					}
					if(length(hypotheses$llist)==0) {
						hypotheses$llist = list(event.map$events.name);
					}
					else {
						hypotheses$llist = list(c(unlist(hypotheses$llist),event.map$events.name));
					}
					for (j in 1:length(event.map$events.name)) {
						function.name = paste(function.name,"_",event.map$events.name[j],sep="");
						function.inputs = c(function.inputs,event.map$events.name[j]);
					}
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
				function.inputs = c(function.inputs,clauses[[i]]$function.name);
			}
		}
		result = list(curr_dataset=curr_dataset,hypotheses=hypotheses,function.name=function.name,function.inputs=function.inputs);
		#save the new edges
		for(k in 1:length(result$function.inputs)) {
			lifting.edges = rbind(lifting.edges,c(result$function.inputs[[k]],result$function.name));
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
	#look for the global variables named lifting.dataset and lifting.annotations
	dataset = lifting.dataset;
	annotations = lifting.annotations;
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations,"AND",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
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
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.inputs=function.inputs);
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
	#look for the global variables named lifting.dataset and lifting.annotations
	dataset = lifting.dataset;
	annotations = lifting.annotations;
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations,"OR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		#evaluate the OR operator
		formula = rep(0,nrow(dataset));
		for (i in 1:nrow(dataset)) {
			formula[i] = sum(curr_dataset[i,]);
			if(formula[i]>1) {
				formula[i] = 1;
			}
		}
		formula = as.integer(formula);
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.inputs=function.inputs);
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
	#look for the global variables named lifting.dataset and lifting.annotations
	dataset = lifting.dataset;
	annotations = lifting.annotations;
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations,"XOR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		#evaluate the XOR operator
		formula = rep(0,nrow(dataset));
		for (i in 1:nrow(dataset)) {
			formula[i] = sum(curr_dataset[i,]);
			if(formula[i]>1) {
				formula[i] = 0;
			}
		}
		formula = as.integer(formula);
		result = list(formula=formula,hypotheses=hypotheses,function.name=function.name,function.inputs=function.inputs);
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA);
}

#### end of file -- lifting.operators.R
