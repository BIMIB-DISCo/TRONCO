#### lifting.operators.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Utility function to add the hypotheses
aux.log = function( dataset, annotations, function.name, ... ) {
  
	if(!is.null(dataset) && !is.null(annotations) && length(list(...)) > 0) {
		clauses = list(...)
		curr_dataset = array(0, c(nrow(dataset), length(clauses)))
		hypotheses = list()
		function.inputs = list()
    
    ### test
    #print('aux log clauses')
    #print(clauses)
    
		for (i in 1:length(clauses)) {
			# if the clause is given by name, get the column from the dataset..
      
      ### test
      #print('typeof cl 1')
      #print(typeof(clauses[[1]]))
      
      
			if(typeof(clauses[[i]]) == "character") {
				
			  col.num = -1
        # if I have the label, get the column in the dataset for this event
				if(length(clauses[[i]]) == 1) {
          ### print
          #print('arg of emap with l = 1')
          #print(c(clauses[[i]],"*"))
          
          
					event.map = emap(c(clauses[[i]],"*"), dataset, annotations)
					col.num = event.map$col.num
				} else if(length(clauses[[i]]) == 2) {
          
				  #print('arg of emap with l = 2')
				  #print(clauses[[i]])
          
          
					event.map = emap(clauses[[i]], dataset, annotations)
					col.num = event.map$col.num
				} 
        
				if(col.num[1] == -1) {
					stop(paste("[ERR] No events for gene ", paste(clauses[[i]], collapse=', ', sep='')))
				} else {
					curr_dataset[,i] = dataset[,col.num[1]]
					if(length(col.num)>1) {
						curr_dataset = cbind(curr_dataset, dataset[ ,col.num[2:length(col.num)]])
					}
          
					if(length(hypotheses$llist) == 0) {
						hypotheses$llist = list(event.map$events.name)
					} else {
						hypotheses$llist = list(c(unlist(hypotheses$llist), event.map$events.name))
					}
					for (j in 1:length(event.map$events.name)) {
						function.name = paste(function.name,"_",event.map$events.name[j],sep="")
						function.inputs = c(function.inputs,event.map$events.name[j])
					}
				}
			} else {
			  # ..otherwise I already have the column as a vector
				
			  curr_dataset[,i] = clauses[[i]]$formula
        # if it is a list
				if(length(hypotheses$llist) == 0) {
					hypotheses$llist = clauses[[i]]$hypotheses$llist
				} else {
					hypotheses$llist = list(c(unlist(clauses[[i]]$hypotheses$llist), unlist(hypotheses$llist)))
				}
				function.name = paste(function.name, "_", clauses[[i]]$function.name, sep="")
				function.inputs = c(function.inputs, clauses[[i]]$function.name)
			}
		}
    
    		# print(curr_dataset)
		result = list(curr_dataset=curr_dataset,
                  hypotheses=hypotheses,
                  function.name=function.name,
                  function.inputs=function.inputs,
                  tests = pairwise.fisher.test(curr_dataset))
    
		#save the new edges
		for(k in 1:length(result$function.inputs)) {
			lifting.edges = rbind(lifting.edges, c(result$function.inputs[[k]], result$function.name))
			assign("lifting.edges", lifting.edges, envir=.GlobalEnv)
		}
		
		return(result)
	} else {
		stop("[ERR] Either the dataset or the formula not provided! No hypothesis will be created.")
	}
	return(NA)
}

# AND hypothesis
"AND" <-
function( ... ) {
	# look for the global variables named lifting.dataset and lifting.annotations
	dataset = lifting.dataset;
	annotations = lifting.annotations;
	pvalue = lifting.pvalue;
  	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		# get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations, "AND" ,...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		#verify the fisher exact tests
		fisher.pvalues = vector();
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				odds.ratio = curr.test[2]
				p.value = curr.test[1]
	
				# print(p.value)
				# if(p.value > pvalue) stop('[ERR] AND: Fisher pvalue >', pvalue, ' (p=', p.value, ') - statistics for 2 genes is not significant.') 			
				# if(odds.ratio <= 0) stop('[ERR] AND: Odds ratio <= 0 (', odds.ratio, ') - 2 of these genes suggest an exclusivity trend.')
	
				# if(curr.test[2]<=0) {
					# curr.pvalue = 1;
				# }
				# else {
					# curr.pvalue = curr.test[1];
				# }
				# if(curr.pvalue>pvalue) {
					# stop(paste("[ERR] Found an invalid pattern with pvalue: ", toString(curr.pvalue), sep=''))
				# }
				fisher.pvalues = append(fisher.pvalues, p.value)
			}
		}
		#print(fisher.pvalues)
		
		# evaluate the AND operator
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
		result = list(formula=formula,hypotheses=hypotheses, function.name=function.name, function.inputs=function.inputs);
		return(result);
	}
	else {
		stop("Either the dataset or the formula not provided! No hypothesis will be created.");
	}
	return(NA)
}

# OR hypothesis
"OR" <-
function( ... ) {
	#look for the global variables named lifting.dataset and lifting.annotations
	dataset = lifting.dataset;
	annotations = lifting.annotations;
	pvalue = lifting.pvalue;
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations,"OR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		#verify the fisher exact tests
		fisher.pvalues = vector();
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				p.value = 1 - curr.test[1];
	
				# print(curr.test)
				# if(p.value > pvalue) stop('[ERR] OR: Fisher pvalue >', pvalue, ' (p=', curr.pvalue, ') - statistics for 2 genes is not significant.') 				
				# if(curr.pvalue>pvalue) {
					# stop(paste("[ERR] Found an invalid pattern with pvalue: ", toString(curr.pvalue), sep=''))
				# }
				fisher.pvalues = append(fisher.pvalues, p.value)
			}
		}
		#print(fisher.pvalues)
		
		# evaluate the OR operator
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
	pvalue = lifting.pvalue;
	if(!is.null(dataset) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the formula from the dataset
		result = aux.log(dataset,annotations,"XOR",...);
		curr_dataset = result$curr_dataset;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		#verify the fisher exact tests
		fisher.pvalues = vector();
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				odds.ratio = curr.test[2]
				p.value = curr.test[1]
		
				# if(p.value > pvalue) stop('[ERR] XOR: Fisher pvalue >', pvalue, '(p=', p.value, ') - statistics for 2 genes is not significant.') 			
				# if(odds.ratio >= 0) stop('[ERR] XOR: Odds ratio >= 0 (', odds.ratio, ') - 2 of these genes suggest a co-occurrence trend.')


# # 				if(curr.test[2]>=0) {
					# curr.pvalue = 1;
				# }
				# else {
					# curr.pvalue = curr.test[1];
				# }
				# if(curr.pvalue>pvalue) {
					# stop(paste("[ERR] Found an invalid pattern with pvalue: ", toString(curr.pvalue), sep=''))
				# }
				fisher.pvalues = append(fisher.pvalues, p.value)
			}
		}
		#print(fisher.pvalues)
		
		# evaluate the XOR operator
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

testing = function(data, g1, g2) {

	# Dataframe di tutto il dataset
	df = data.frame(row.names=as.samples(data))
	df$x = rowSums(as.gene(data, genes=g1))
	df$y = rowSums(as.gene(data, genes=g2))
	
	# Lifting xor
	df$xor = df$x + df$y
	df$xor[ df$xor > 1] = 0
	
	# Lifting or
	df$or = df$x + df$y
	df$or[ df$or > 1] = 1
	
	# Lifting and
	df$and = df$x + df$y
	df$and[ df$and < 2] = 0
	df$and[ df$and == 2] = 1
	
	# Nomi per accedere successivamente 
	names(df$x) = g1
	names(df$y) = g2
	names(df$xor) = 'xor'
	names(df$or) = 'or'
	names(df$and) = 'and'
	
	cat('DATASET\n')
	print(df)
	
	# Tabella di contingenza 2x2
	table.xor = rbind(
		c(nrow(df) - sum(df$or), sum(df$or - df$y)),
		c(sum(df$or - df$x), sum(df$and))
	)
	
	colnames(table.xor) = c(paste0('-', g1), paste0('+', g1))
	rownames(table.xor) = c(paste0('-', g2), paste0('+', g2))
	
	cat('\nCATEGORICAL\n')
	print(table.xor)
	
	# Fisher 2-sided
	test = fisher.test(table.xor)
	
	# p-value e log dell’odds ratio
	cat('p-value (2-sided): ', test$p.value, '\n')
	cat('log(odds ratio): ', log(test$estimate['odds ratio']))
}

# performs pairwise exact fisher test
pairwise.fisher.test = function(data) {
	
	# structure to save the results
	results = vector();
	
	if(ncol(data)>1) {
		for(i in 1:ncol(data)) {
			for(j in i:ncol(data)) {
				if(i!=j) {
					
					df = data[,c(i,j)]
					df_x = data[,1]
					df_y = data[,2]
					
					# Lifting xor
					df_xor = df_x + df_y
					df_xor[ df_xor > 1] = 0
					
					# Lifting or
					df_or = df_x + df_y
					df_or[ df_or > 1] = 1
					
					# Lifting and
					df_and = df_x + df_y
					df_and[ df_and < 2] = 0
					df_and[ df_and == 2] = 1
					
					# 2x2 contingency table
					table.xor = rbind(
						c(nrow(df) - sum(df_or), sum(df_or - df_y)),
						c(sum(df_or - df_x), sum(df_and))
					)
					
					# Fisher 2-sided
					test = fisher.test(table.xor)
					
					# save the results
					curr_result = c(test$p.value,log(test$estimate['odds ratio']))
					results = append(results, list(curr_result))
					
					# print('TESTING')
					# print(colnames(data)[i])
					# print(colnames(data)[j])
					# print(df)
					# print(table.xor)
					# print(test)
				}
			}
		}
	}
	
	return(results)
	
}

#### end of file -- lifting.operators.R
