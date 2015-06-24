#### hypothesis.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Add a new hypothesis by creating a new event and adding it to the compliant genotypes
#' @export
"hypothesis.add" <-
  function( data, pattern.label, lifted.pattern, pattern.effect = "*", pattern.cause = "*" ) {

	# save the needed data structures
    if(!is.null(data$genotypes) && !is.null(data$annotations)) {
    		genotypes = data$genotypes;
    		annotations = data$annotations;
    }
    else {
    		genotypes = NULL;
    		annotations = NULL;
    }
    if(!is.null(data$hypotheses)) {
    		hypotheses = data$hypotheses;
    }
    else {
    		hypotheses = NA;
    }
    
    # add the hypothesis only if all the inputs are correctly provided
    if(!is.null(genotypes) && !is.null(annotations)) {
    		
    		# the Boolean functions look for a set of global variables
      	# if there are already global variables named as the ones used here, make the backup of them
      	do.roll.back.lifting.genotypes = FALSE;
      	do.roll.back.lifting.annotations = FALSE;
      	do.roll.back.lifting.edges = FALSE;
      	do.roll.back.fisher.pvalues = FALSE;
      
      	# I need a global variable to save the genotypes of the lifted pattern
      	# if there is already a global variable named lifting.genotypes, make the backup of it
      	if(exists("lifting.genotypes")) {
      		roll.back.lifting.genotypes = lifting.genotypes;
        		do.roll.back.lifting.genotypes = TRUE;
      	}
      	assign("lifting.genotypes",genotypes,envir=.GlobalEnv);
      	
      	# I need a global variable to save the annotations of the lifted pattern
	    # if there is already a global variable named lifting.annotations, make the backup of it
	    if(exists("lifting.annotations")) {
	    		roll.back.lifting.annotations = lifting.annotations;
	      	do.roll.back.lifting.annotations = TRUE;
	    }
	    assign("lifting.annotations",annotations,envir=.GlobalEnv);
	    
	    # I need a global variable to save the edges of the lifted pattern
      	# if there is already a global variable named lifting.edges, make the backup of it
      	if(exists("lifting.edges")) {
      		roll.back.lifting.edges = lifting.edges;
      		do.roll.back.lifting.edges = TRUE;
      	}
      	assign("lifting.edges",NULL,envir=.GlobalEnv);
	    
	    # I need a global variable to save the pvalues of the lifted pattern
      	# if there is already a global variable named fisher.pvalues, make the backup of it
      	if(exists("fisher.pvalues")) {
      		roll.back.fisher.pvalues = fisher.pvalues;
      		do.roll.back.fisher.pvalues = TRUE;
      	}
      	assign("fisher.pvalues",NULL,envir=.GlobalEnv);
      	
      	# save the lifted genotypes and its hypotheses for the current pattern
      	curr_pattern = lifted.pattern$pattern;
      	curr_hypotheses = lifted.pattern$hypotheses;
      	curr_pvalues = fisher.pvalues;
      	
      	# save the edges of the lifted pattern
      	hstructure = lifting.edges;
      	
      	# roll back to the previous value of the global variable lifting.genotypes if any or remove it
      	if(do.roll.back.lifting.genotypes) {
      		assign("lifting.genotypes",roll.back.lifting.genotypes,envir=.GlobalEnv);
      	}
      	else {
      		rm(lifting.genotypes,pos=".GlobalEnv");
      	}
      	
      	# roll back to the previous value of the global variable lifting.annotations if any or remove it
      	if(do.roll.back.lifting.annotations) {
      		assign("lifting.annotations",roll.back.lifting.annotations,envir=.GlobalEnv);
      	}
      	else {
      		rm(lifting.annotations,pos=".GlobalEnv");
      	}
      	
      	# roll back to the previous value of the global variable lifting.edges if any or remove it
      	if(do.roll.back.lifting.edges) {
      		assign("lifting.edges",roll.back.lifting.edges,envir=.GlobalEnv);
      	}
      	else {
      		rm(lifting.edges,pos=".GlobalEnv");
      	}
      	
      	# roll back to the previous value of the global variable fisher.pvalues if any or remove it
      	if(do.roll.back.fisher.pvalues) {
      		assign("fisher.pvalues",roll.back.fisher.pvalues,envir=.GlobalEnv);
      	}
      	else {
      		rm(fisher.pvalues,pos=".GlobalEnv");
      	}
		
      	# set the hypotheses number
      	if(!is.na(hypotheses[1])) {
      		num.hypotheses = hypotheses$num.hypotheses;
      	}
      	else {
      		num.hypotheses = 0;
      	}
      	
      	# * is a special pattern.effect which indicates to use all the events as effects for this pattern
      	is.to.all.effects = FALSE;
      	if(pattern.effect[[1]][1]=="*") {
      		
      		pattern.effect = colnames(genotypes)[1:(length(colnames(genotypes))-num.hypotheses)];
      		
      		# any event can not be both causes and effects for the pattern to be well-formed
      		pattern.effect = list(pattern.effect[-which((pattern.effect%in%unlist(curr_hypotheses$llist)))]);
      		is.to.all.effects = TRUE;
			
      		if(length(pattern.effect)==0) {
      			stop(paste("[ERR] Missing list of effects to test or wildcard \'*\'.", sep=''));
      		}
      		
      	}
      	
      	# check the pattern to be well-formed
      	all.col.nums = vector();
      	if(length(pattern.effect)==0) {
      		stop(paste("[ERR] Missing list of effects or wildcard \'*\'.", sep=''));
      	}
      	else {
      		# check the effects of the pattern to be well-formed
      		for (i in 1:length(pattern.effect)) {
      			
      			curr.pattern.effect = pattern.effect[[i]];
      			if(is.to.all.effects==FALSE) {
      				col.num = -1;
      				if(length(curr.pattern.effect)==1) {
      					event.map = emap(c(curr.pattern.effect,"*"),genotypes,annotations);
      					col.num = event.map$col.num;
      					events.name = event.map$events.name;
      				}
      				else if(length(curr.pattern.effect)==2) {
      					event.map = emap(curr.pattern.effect,genotypes,annotations);
      					col.num = event.map$col.num;
      					events.name = event.map$events.name;
      				}
      			}
      			else {
      				col.num = which(colnames(genotypes)%in%curr.pattern.effect);
      				if(length(col.num)==0) {
      					col.num = -1;
      				}
      				events.name = curr.pattern.effect;
      			}
      			
      			# check the effect to be a valid event
      			if(col.num[1]==-1) {
      				stop(paste("[ERR] Unknown gene among effects: \"", curr.pattern.effect, "\".",sep=''));
      			}
      			all.col.nums = append(all.col.nums,col.num);
      			
      			# check the pattern to be well-formed
      			# if the effect is in the pattern, the pattern is not well-formed
      			if(length(which(unlist(curr_hypotheses$llist)%in%events.name))>0) {
      				stop(paste("[ERR] Bad formed pattern, event \"", curr.pattern.effect, "\" yields a loop.",sep=''));
      			}
      			
      		}
      	}
      	
      	# look for duplicated effects in the pattern
      	if(anyDuplicated(all.col.nums)>0) {
      		stop(paste("[ERR] Bad formed pattern, duplicated events ", 
                   paste(pattern.effect[duplicated(pattern.effect)], collapse=', ', sep=''),
                   "within effects.", sep=''));
        }
        
        # check that the we are not duplicating any name by adding the new pattern
        if(length(which(colnames(genotypes)==pattern.label))>0) {
        		stop(paste("[ERR] This pattern already exists.", sep=''));
        	}
		
        	# add the pattern to the genotypes
        	genotypes = cbind(genotypes,curr_pattern);
        	
        	# check that the pattern is valid according to Suppes' theory
        	# compute the observed and observed joint probabilities
        	pair.count <- array(0, dim=c(ncol(genotypes),ncol(genotypes)));
        	# compute the probabilities on the genotypes
        	for(i in 1:ncol(genotypes)) {
        		for(j in 1:ncol(genotypes)) {
        			val1 = genotypes[ ,i];
        			val2 = genotypes[ ,j];
        			pair.count[i,j] = (t(val1)%*%val2);
        		}
        	}
        	# marginal.probs is an array of the observed marginal probabilities
        	marginal.probs <- array(as.matrix(diag(pair.count)/nrow(genotypes)),dim=c(ncol(genotypes),1));
        	# joint.probs is an array of the observed joint probabilities
        	joint.probs <- as.matrix(pair.count/nrow(genotypes));
        	# check that the probability of the pattern is in (0,1)
        	if(marginal.probs[ncol(genotypes)]==0 || marginal.probs[ncol(genotypes)]==1) {
        		stop(paste("[ERR] The pattern has marginal probability ", marginal.probs[ncol(genotypes)], 
                   ", which should be in (0,1).", sep=''));
        }
        
        # check that the pattern does not duplicate any existing column
        i = ncol(genotypes);
        for(j in 1:ncol(genotypes)) {
        		# if the edge is valid, i.e., not self cause
        		if(i!=j) {
        			#if the two considered events are not distinguishable
        			if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) {
        				stop(paste("[ERR] Pattern duplicates ", paste(as.events(data)[j, ], collapse=' ', sep=''), 
                       ".", sep=''));
                }
            }
        }
        
        # * is a special pattern.cause which indicates to use all the events as causes for this pattern
      	is.to.all.causes = FALSE;
      	if(pattern.cause[[1]][1]=="*") {
      		
      		pattern.cause = colnames(genotypes)[1:(length(colnames(genotypes))-num.hypotheses-1)];
      		
      		# any event can not be both causes and effects for the pattern to be well-formed
      		pattern.cause = list(pattern.cause[-which((pattern.cause%in%unlist(curr_hypotheses$llist)))]);
      		is.to.all.causes = TRUE;
      		
      	}
      	
      	# check the pattern to be well-formed
      	all.col.nums = vector();
      	if(length(pattern.cause)>0) {
      		
      		# check the causes of the pattern to be well-formed
      		for (i in 1:length(pattern.cause)) {
      			
      			curr.pattern.cause = pattern.cause[[i]];
      			if(is.to.all.causes==FALSE) {
      				col.num = -1;
      				if(length(curr.pattern.cause)==1) {
      					event.map = emap(c(curr.pattern.cause,"*"),genotypes,annotations);
      					col.num = event.map$col.num;
      					events.name = event.map$events.name;
      				}
      				else if(length(curr.pattern.cause)==2) {
      					event.map = emap(curr.pattern.cause,genotypes,annotations);
      					col.num = event.map$col.num;
      					events.name = event.map$events.name;
      				}
      			}
      			else {
      				col.num = which(colnames(genotypes)%in%curr.pattern.cause);
      				if(length(col.num)==0) {
      					col.num = -1;
      				}
      				events.name = curr.pattern.cause;
      			}
      			
      			# check the cause to be a valid event
      			if(col.num[1]==-1) {
      				stop(paste("[ERR] Unknown gene among causes: \"", curr.pattern.cause, "\".",sep=''));
      			}
      			all.col.nums = append(all.col.nums,col.num);
      			
      			# check the pattern to be well-formed
      			# if the cause is in the pattern, the pattern is not well-formed
      			if(length(which(unlist(curr_hypotheses$llist)%in%events.name))>0) {
      				stop(paste("[ERR] Bad formed pattern, event \"", curr.pattern.cause, "\" yields a loop.",sep=''));
      			}
      			
      		}
      	}
      	
      	# look for duplicated causes in the pattern
      	if(anyDuplicated(all.col.nums)>0) {
      		stop(paste("[ERR] Bad formed pattern, duplicated events ", 
                   paste(pattern.cause[duplicated(pattern.cause)], collapse=', ', sep=''),
                   "within causes.", sep=''));
        }
        
        # now I can finally add the hypothesis
        colnames(genotypes)[ncol(genotypes)] = pattern.label;
        if(is.na(hypotheses[1])) {
        		hypotheses = list();
        	}
        	hypotheses$num.hypotheses = num.hypotheses + 1;
		
        # add the new hypothesis in the annotations
        annotations = rbind(data$annotations,c("Pattern", pattern.label));
        rownames(annotations)[nrow(annotations)] = pattern.label;
        
        # add the color of the type "Hypothesis" is not already defined
        if(any(rownames(data$types)=="Pattern")==FALSE) {
			types = rbind(data$types, 'slateblue');
			rownames(types)[nrow(types)] = "Pattern";
			data$types = types;
        }
        	
        	# create the list of added hypotheses
        	if(length(hypotheses$hlist)==0) {
        		hypotheses$hlist = vector();
        	}
        	# add the new hypothesis to the list
        	for (i in 1:length(pattern.effect)) {
        		curr.pattern.effect = pattern.effect[[i]];
        		if(is.to.all.effects==FALSE) {
        			if(length(curr.pattern.effect)==1) {
        				event.map = emap(c(curr.pattern.effect,"*"),genotypes,annotations);
        				col.num = event.map$col.num;
        			}
        			else if(length(curr.pattern.effect)==2) {
        				event.map = emap(curr.pattern.effect,genotypes,annotations);
        				col.num = event.map$col.num;
        			}
        		}
        		else {
        			col.num = which(colnames(genotypes)%in%curr.pattern.effect);
        			if(length(col.num)==0) {
        				col.num = -1;
        			}
        		}
        		for (j in 1:length(col.num)) {
        			hypotheses$hlist = rbind(hypotheses$hlist,t(c(colnames(genotypes)[ncol(genotypes)],colnames(genotypes)[col.num[j]])));
        		}
        		if(is.null(colnames(hypotheses$hlist))) {
        			colnames(hypotheses$hlist) = c("cause","effect");
        		}
        	}
        	
        # add the causes of the new hypothesis to the list
        if(length(pattern.cause)>0) {
        		for (i in 1:length(pattern.cause)) {
        			curr.pattern.cause = pattern.cause[[i]];
        			if(is.to.all.causes==FALSE) {
        				if(length(curr.pattern.cause)==1) {
        					event.map = emap(c(curr.pattern.cause,"*"),genotypes,annotations);
        					col.num = event.map$col.num;
        				}
        				else if(length(curr.pattern.cause)==2) {
        					event.map = emap(curr.pattern.cause,genotypes,annotations);
        					col.num = event.map$col.num;
        				}
        			}
        			else {
        				col.num = which(colnames(genotypes)%in%curr.pattern.cause);
        				if(length(col.num)==0) {
        					col.num = -1;
        				}
        			}
        			for (j in 1:length(col.num)) {
        				hypotheses$hlist = rbind(hypotheses$hlist,t(c(colnames(genotypes)[col.num[j]],colnames(genotypes)[ncol(genotypes)])));
        			}
        			if(is.null(colnames(hypotheses$hlist))) {
        				colnames(hypotheses$hlist) = c("cause","effect");
        			}
        		}
        }
        	
        	# create the list of hypotheses' structures
        	if(length(hypotheses$hstructure)==0) {
        		hypotheses$hstructure = new.env(hash=TRUE,parent=emptyenv());
        	}
        	hypotheses$hstructure[[pattern.label]] = get.lifted.pattern(hstructure);
        	
        	# add the atoms of the hypothesis
        if(length(hypotheses$patterns)==0) {
			hypotheses$patterns = list();
        }
        hypotheses$patterns[pattern.label] = lifted.pattern$hypotheses$llist;
        
        #add the hypotheses of the atoms
        if(length(hypotheses$atoms)==0) {
        		hypotheses$atoms = vector(mode="list",length=(ncol(genotypes)-hypotheses$num.hypotheses));
        		names(hypotheses$atoms) = colnames(genotypes)[1:(ncol(genotypes)-hypotheses$num.hypotheses)];
        	}
        	atoms.in.pattern = which(names(hypotheses$atoms)%in%unlist(hypotheses$patterns[pattern.label]));
        	if(length(atoms.in.pattern)>0) {
        		for (i in 1:length(atoms.in.pattern)) {
        			hypotheses$atoms[[atoms.in.pattern[i]]] = append(hypotheses$atoms[[atoms.in.pattern[i]]], pattern.label);
        		}
        	}
        
        #add the fisher pvalues
        if(length(hypotheses$pvalues)==0) {
        		hypotheses$pvalues = vector();
        }
        hypotheses$pvalues = append(hypotheses$pvalues,list(curr_pvalues))
        names(hypotheses$pvalues)[length(hypotheses$pvalues)] = pattern.label;
        	
        	#return the new (compliant) data structure as result
        	data$genotypes = genotypes;
        data$annotations = annotations;
        	data$hypotheses = hypotheses;
        	
        	return(data);
    }
    else {
    		stop("[ERR] Missing genotypes or pattern.");
    		
    	}
    	
    return(NA);
    
}

#### end of file -- hypothesis.add.R

# resolve the ellipsis for the effects
hypothesis.lifted.effects = function( ... ) {
	return(list(...));
}

#' @export hypothesis.add.group
hypothesis.add.group = function(x, 
                                FUN, 
                                group, 
                                pattern.cause = '*',
                                pattern.effect = '*',
                                dim.min = 2,
                                dim.max = length(group),
                                min.prob = 0) 
{
	op = deparse(substitute(FUN))

	#print(length(unlist(group)))
	#print(unlist(group))

	#effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste0("c('", paste(pattern.effect, collapse = "', '"), "')")
  cause = paste0("c('", paste(pattern.cause, collapse = "', '"), "')")

  # print(effect)
  # print(cause)
	
	ngroup = length(group)
	if (ngroup < 2) {
		warning("No hypothesis will be created for groups with less than 2 elements.")
		return(x)
	}

	cat("*** Adding Group Hypotheses\n")
	cat(' Group:', paste(group, collapse = ", ", sep = ""), '\n')
	cat(' Function:', op, '\n')
   cat(' Cause:', paste(pattern.cause, collapse=", "), '; ')
   cat(' Effect:', paste(pattern.effect, collapse=", "), '.\n')
	flush.console()

	if(min.prob > 0)
	{
		cat('\nFiltering genes within the group with alteration frequency below', min.prob, '\n')
		
		temp = events.selection(x, filter.in.names = group)
		
		#show(temp)
		temp = as.alterations(temp)
		temp = events.selection(temp, filter.freq = min.prob)
		
		group = as.genes(temp)
		cat('New group:', paste(group, collapse = ", ", sep = ""), '\n')
	}
	
	ngroup = length(group)
	if (ngroup < 2) {
		warning("No hypothesis will be created for groups with less than 2 elements.")
		return(x)
	}
	
	hom.group = lapply(group, function(g, x) {
		if (nevents(x, genes = g) > 1) 
			T
		else F
	}, x)
	hom.group = group[unlist(hom.group)]

	gene.hom = function(g, h) {
		if (g %in% h) 
		{
			if( any(rowSums(as.gene(x, genes = g)) > 1) ) return(paste0("OR('", g, "')"))
			else return(paste0("XOR('", g, "')"))
		}
		return(paste0("'", g, "'"))
	}

	max.groupsize = min(dim.max, ngroup)
	min.groupsize = max(2, dim.min)
	if(dim.min > dim.max) stop('ERROR - dim.min > dim.max')
	if(min.groupsize > max.groupsize) stop('ERROR - min.groupsize > max.groupsize')
	
	if (length(hom.group) > 0) 
		cat("Genes with multiple events: ", paste(unlist(hom.group), collapse=', ', sep=''), "\n")
	
	error.summary = data.frame()

	# Get an analytical pattern... !
	tot.patterns = 0
	for (i in min.groupsize:max.groupsize) tot.patterns = tot.patterns + ncol(combn(unlist(group), i))
	
	# create a progress bar
	cat('Generating ', tot.patterns ,'patterns [size: min =', max.groupsize,' -  max =', max.groupsize, '].\n')
		
	# pb <- txtProgressBar(0, tot.patterns, style = 3)
	flush.console()

	pbPos = 0
	for (i in min.groupsize:max.groupsize) {
		gr = combn(unlist(group), i)
	
		#print(gr)
		
		
		for (j in 1:ncol(gr)) {
			genes = as.list(gr[, j])

			#start the progress bar
			pbPos = pbPos + 1
			# setTxtProgressBar(pb, pbPos)

			hypo.name = paste(unlist(genes), sep = "_", collapse = "_")
			hypo.genes = paste(lapply(genes, function(g, hom.group) {
				gene.hom(g, hom.group)
			}, hom.group), collapse = ", ")

			# print(hypo.genes)
			# print(hom.group)

				hypo.add = paste0("hypothesis.add(x, ", 
                          "pattern.label= '", op, "_", hypo.name, "', ",
                          "lifted.pattern= ", op, "(", hypo.genes, "), ",
                          "pattern.effect=", effect, ", ",
                          "pattern.cause=", cause, ")")

				# cat('*** Evaluating ', hypo.add, '\n')
				
				err = tryCatch({
					x = eval(parse(text = hypo.add))
				}, error = function(cond) {
					# print(cond)
					m = paste("Error on", hypo.add, ".\n", cond)
					code = strsplit(as.character(cond), " ")[[1]]
					idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

					return(
						data.frame(
							pattern = paste(unlist(genes), collapse = ", ", sep = ""), 
							error = paste(code[idx.errcode:length(code)], collapse = " ")
							))

				}, warning = function(cond) {
					m = paste("Warning on", hypo.add, ".\n", cond)
					return(genes)
				})
				# Dummy errors detection
				if (!("genotypes" %in% names(err))) 
					error.summary = rbind(error.summary, err)
			}
	}
	# close progress bar
	# close(pb)

	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " genes pattern could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible patterns.\n")

	return(x)
}

#' @export
hypothesis.add.homologous = function(x, 
                                     pattern.cause = '*',
                                     pattern.effect = '*',
                                     genes = as.genes(x),
                                     FUN = OR) 
{
	# in questa funzione, per ogni gene che ha più di un tipo di alterazione
	# aggiungo l'OR

	#effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	#effect = paste(effect, collapse = ", ")
	op = deparse(substitute(FUN))

	hom.group = lapply(genes, function(g, x) {
		if (nevents(x, genes = g) > 1) 
			T
		else F
	}, x)
	hom.group = genes[unlist(hom.group)]

	if (length(hom.group) == 0) {
		warning("No genes with multiple events.")
		return(x)
	}

	cat("*** Adding hypotheses for Homologous Patterns\n")
	cat(' Genes:', paste(hom.group, collapse = ", ", sep = ""), '\n')
	cat(' Function:', op, '\n')
    cat(' Cause:', paste(pattern.cause, collapse=", "), '\n')
	cat(' Effect:', paste(pattern.effect, collapse=", "), '\n')
	flush.console()

  effect = paste0("c('", paste(pattern.effect, collapse = "', '"), "')")
  cause = paste0("c('", paste(pattern.cause, collapse = "', '"), "')")

	if (length(hom.group) == 0) {
		warning("No genes with multiple events.")
		return(x)
	}

	
	# print(length(hom.group))
	# create a progress bar
	pb <- txtProgressBar(0, length(hom.group), style = 3)

	

	error.summary = data.frame()

	for (i in 1:length(hom.group)) {

		#start the progress bar
		setTxtProgressBar(pb, i)

		# Check if the joint probability of homologous events is > 0, if
		# yes the event will be added as 'OR', otherwise 'XOR'
		if( any(rowSums(as.gene(x, genes = hom.group[[i]])) > 1))
		  FUN = 'OR'
		else
      FUN = 'XOR'				

		hypo.add = paste0("hypothesis.add(x, ",
                      "pattern.label= '", FUN, "_", hom.group[[i]], "', ",
                      "lifted.pattern= ", FUN, "('", hom.group[[i]], "'), ",
                      "pattern.cause= ",  cause, ", ",
                      "pattern.effect=", effect, ")")

    # cat('*** Evaluating ', hypo.add, '\n')

		err = tryCatch({
			x = eval(parse(text = hypo.add))
		}, error = function(cond) {
			m = paste("Error on", hypo.add, ".\n", cond)
			code = strsplit(as.character(cond), " ")[[1]]
			idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

			return(data.frame(pattern = paste(unlist(hom.group[[i]]), collapse = ", ", sep = ""), error = paste(code[idx.errcode:length(code)], collapse = " ")))

		}, warning = function(cond) {
			m = paste("Warning on", hypo.add, ".\n", cond)
			return(genes)
		})
		# Dummy errors detection
		if (!("genotypes" %in% names(err))) 
			error.summary = rbind(error.summary, err)
	}

	# close progress bar
	close(pb)

	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " patterns could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible gene patterns.\n")


	return(x)

}

#' @import igraph
hypotheses.expansion <- function(input_matrix, 
                                 map = list(),
                                 hidden_and = T,
                                 expand = T,
                                 events = NULL,
                                 # conf_matrix = NULL,
                                 skip.disconnected = TRUE
                                 ) {
  

  suppressMessages(library(igraph))
  
  
  # get node list
  node_list <- colnames(input_matrix)
  #print('input matrix')
  #print(input_matrix)

  #print('map')
  #print(names(map))
  
  # cut input matrix
  num_hypos = 0
  if(length(map) > 0) {
    num_hypos = Reduce(sum, lapply(ls(map), function(x, y){if(x %in% y)return(1)}, y=node_list))
  }
 #print('num_hypos')
  #print(num_hypos)

  margin = length(node_list) - num_hypos
  hypos_new_name = list()
  
  # check if there are hypotheses
  if (num_hypos == 0 || !expand) {
    # if no hypos do nothings..
    min_matrix = input_matrix
  } else {
    # ..else expand them

    #print('input matrix')
    #print(input_matrix)

    min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]

    #print('min matrix')
    #print(min_matrix)
    
    # create graph from matrix
    min_graph = graph.adjacency(min_matrix)
    
    
    # foreach hypothesis
    # print(ls(map))
    # for (h in ls(map)) {
    # print(h)
    # }
      
    for (h in ls(map)) {
      #print('hypo to expand')
      #print(h)
      if (! h %in% node_list) {
        next
      }
      
      # eros! please give me the transposed matrix
      hypo = map[[h]]
      
      # create graph from hypo
      hypo_graph = graph.adjacency(hypo)

      # name of this node
      h_mat <- rowSums(get.adjacency(hypo_graph, sparse=FALSE))
      
      initial_node <- names(h_mat)[which(h_mat==0)]
      hypos_new_name[initial_node] = h

      #print(paste("new name:", initial_node))

      # change names in confidence matrix according to hypotesis
      # if(!is.null(conf_matrix)) {
      #   rownames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
      #   colnames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
      # }

      display.up = FALSE
      if (length(which(input_matrix[, h] == 1)) != 0) {
        display.up = TRUE
      }

      display.down = FALSE
      if (length(which(input_matrix[h, ] == 1)) != 0) {
        display.down = TRUE
      }
      
      # print('***')
      # print(h)
      # print(display.up)
      # print(display.down)

      # display up hypo and reconnect
      if (display.up) {

        hypo_pre = t(hypo)

        node_names = rownames(hypo_pre)
        node_names = lapply(node_names, function(x){ if(is.logic.node(x)) { paste0('UP', x) } else { return(x) }  })



        rownames(hypo_pre) = node_names
        colnames(hypo_pre) = node_names

        # create graph from hypo
        hypo_graph_pre = graph.adjacency(hypo_pre)

        # name of this node
        h_mat_pre <- colSums(get.adjacency(hypo_graph_pre, sparse=FALSE))

        final_node <- names(h_mat_pre)[which(h_mat==0)]
        hypos_new_name[final_node] = h

        # edge to reconstruct
        h_edge <- input_matrix[, h]
        initial_node_up <- names(h_edge)[which(h_edge==1)]

        #print("pre merge")
        #print(get.adjacency(min_graph, sparse=FALSE))

        # add this graph to main graph
        min_graph = graph.union(min_graph, hypo_graph_pre)

        #print("post merge")
        #print(get.adjacency(min_graph, sparse=FALSE))

        # recreate lost edge
        for (node in initial_node_up) {
          #print(paste('new edge:', node, "->", final_node))
          min_graph <- min_graph + edge(node, final_node)
        }

      }
      
      # display down hypo
     # if (display.up || display.down) {
     if (display.down) {

        # edge to reconstruct
        h_edge <- input_matrix[h,]
        final_node <- names(h_edge)[which(h_edge==1)]

        #print("pre merge")
        #print(get.adjacency(min_graph, sparse=FALSE))

        # add this graph to main graph
        min_graph = graph.union(min_graph, hypo_graph)

        #print("post merge")
        #print(get.adjacency(min_graph, sparse=FALSE))
        
      #}

      # reconnect down hypo
      #if (display.down) {
        # print(final_node)
        # recreate lost edge
        for (node in final_node) {
          #print(paste('new edge:', initial_node, "->",node))
          min_graph <- min_graph + edge(initial_node, node)
        }
      }
      
    }
    min_matrix = get.adjacency(min_graph, sparse = F)
    
  }


  # now expand the hidden AND
  #print(min_matrix)
  
  if(hidden_and == F) {
    # sort col and row (igraph wants the same order)
    min_matrix = min_matrix[,order(colnames(min_matrix))]
    min_matrix = min_matrix[order(rownames(min_matrix)),]
    
    # print(min_matrix)
    #if(!is.null(conf_matrix)) {
    #  return(list(min_matrix, hypos_new_name, conf_matrix))
    #}
    return(list(min_matrix, hypos_new_name))
  }
  
  cat('\n*** Expand hidden and:')
  
  and_matrix = NULL
  to_reconnect = list()
  logical_op = list("AND", "OR", "NOT", "XOR", "UPAND", 'UPOR', 'UPXOR')
#   logical_op = list("", "", "", "")
 
  # foreach AND column
  for (col in colnames(min_matrix)) {
    prefix = gsub("_.*$", "", col)
    if ( !(prefix %in% logical_op) && sum(min_matrix[,col]) > 1 ) {
      # not logical operator and colsum > 1 there is a hidden AND and something has to be done..
      
      # remember to reconnect the fake and to this node
      to_reconnect = append(to_reconnect, col)
      # append a colum from the old matrix..
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      # and give her a new name based on the old one
      new_col_name = paste("*", col, sep="_")
      colnames(and_matrix)[pos] = new_col_name
      
      # append a 0 columl to the matrix..
      and_matrix = cbind(and_matrix, matrix(0,nrow = nrow(and_matrix), ncol = 1))
      pos = ncol(and_matrix)
      # and give her the old name
      colnames(and_matrix)[pos] = col
      
      # now do the same to conf_matrix
      if(!is.null(conf_matrix)) {
        conf_matrix = cbind(conf_matrix, conf_matrix[, col])
        pos = ncol(conf_matrix)
        colnames(conf_matrix)[pos] = new_col_name
      }
      
    } else {
      # ..else add the row taken from the old matrix and set the correct colname
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      colnames(and_matrix)[pos] = col
    }
  }
  
  # now reconnect AND node to his gene (AND_Gene7 -> Gene7)
  for(row in to_reconnect) {
    and_matrix = rbind(and_matrix, matrix(0, ncol=ncol(and_matrix), nrow = 1))
    pos = nrow(and_matrix)
    rownames(and_matrix)[pos] = paste0("*_", row)
    and_matrix[paste0("*_", row),row] = 1
  }
  
  # cat(' done')
  
  # sort col and row (igraph wants the same order)
  and_matrix = and_matrix[,order(colnames(and_matrix))]
  and_matrix = and_matrix[order(rownames(and_matrix)),]
  
  # print(and_matrix)
  if(!is.null(conf_matrix)) {
    return(list(and_matrix, hypos_new_name, conf_matrix))
  }
  return(list(and_matrix, hypos_new_name))
}

#### lifting.operators.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Utility function to add the hypotheses
aux.log = function( genotypes, annotations, function.name, ... ) {
  
	if(!is.null(genotypes) && !is.null(annotations) && length(list(...)) > 0) {
		
		clauses = list(...)
		curr_genotypes = array(0, c(nrow(genotypes), length(clauses)))
		hypotheses = list()
		function.inputs = list()
		fisher.tests = vector()
    
		for (i in 1:length(clauses)) {
			
			# if the clause is given by name, get the column from the genotypes
			if(typeof(clauses[[i]]) == "character") {
				
			  	col.num = -1
			  	# if I have the label, get the column in the genotypes for this event
				if(length(clauses[[i]]) == 1) {
					event.map = emap(c(clauses[[i]],"*"), genotypes, annotations)
					col.num = event.map$col.num
				}
				else if(length(clauses[[i]]) == 2) {
					event.map = emap(clauses[[i]], genotypes, annotations)
					col.num = event.map$col.num
				}
        
				if(col.num[1] == -1) {
					stop(paste("[ERR] No events for gene ", paste(clauses[[i]], collapse=', ', sep='')))
				}
				else {
					curr_genotypes[,i] = genotypes[,col.num[1]]
					if(length(col.num)>1) {
						curr_genotypes = cbind(curr_genotypes, genotypes[ ,col.num[2:length(col.num)]])
					}
          
					if(length(hypotheses$llist) == 0) {
						hypotheses$llist = list(event.map$events.name)
					}
					else {
						hypotheses$llist = list(c(unlist(hypotheses$llist), event.map$events.name))
					}
					for (j in 1:length(event.map$events.name)) {
						function.name = paste(function.name,"_",event.map$events.name[j],sep="")
						function.inputs = c(function.inputs,event.map$events.name[j])
					}
				}
				
			}
			else {
				# otherwise I already have the column as a vector
				curr_genotypes[,i] = clauses[[i]]$pattern
				# if it is a list
				if(length(hypotheses$llist) == 0) {
					hypotheses$llist = clauses[[i]]$hypotheses$llist
				}
				else {
					hypotheses$llist = list(c(unlist(clauses[[i]]$hypotheses$llist),unlist(hypotheses$llist)))
				}
				function.name = paste(function.name, "_", clauses[[i]]$function.name, sep="")
				function.inputs = c(function.inputs, clauses[[i]]$function.name)
				
			}
		}
		
		result = list(curr_genotypes=curr_genotypes,
                  hypotheses=hypotheses,
                  function.name=function.name,
                  function.inputs=function.inputs,
                  tests = pairwise.fisher.test(curr_genotypes))
    
		# save the new edges
		for(k in 1:length(result$function.inputs)) {
			lifting.edges = rbind(lifting.edges, c(result$function.inputs[[k]], result$function.name))
			assign("lifting.edges", lifting.edges, envir=.GlobalEnv)
		}
		
		return(result)
		
	}
	else {
		stop("[ERR] Either the genotypes or the pattern not provided! No hypothesis will be created.")
	}
	return(NA)
}

# AND hypothesis
#' @export
"AND" <-
function( ... ) {
	# look for the global variables named lifting.genotypes and lifting.annotations
	genotypes = lifting.genotypes;
	annotations = lifting.annotations;
  	if(!is.null(genotypes) && !is.null(annotations) && length(list(...))>0) {
		# get the vector of the clauses of the pattern from the genotypes
		result = aux.log(genotypes,annotations, "AND" ,...);
		curr_genotypes = result$curr_genotypes;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		#save the fisher exact tests
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				odds.ratio = curr.test[2]
	
				if(curr.test[2]<=0) {
					curr.pvalue = 1;
				}
				else {
					curr.pvalue = curr.test[1];
				}
				
				fisher.pvalues = append(fisher.pvalues,curr.pvalue)
			}
			assign("fisher.pvalues", fisher.pvalues, envir=.GlobalEnv)
		}
		
		# evaluate the AND operator
		pattern = rep(0,nrow(genotypes));
		for (i in 1:nrow(genotypes)) {
			pattern[i] = sum(curr_genotypes[i,]);
			if(pattern[i]<ncol(curr_genotypes)) {
				pattern[i] = 0;
			}
			else {
				pattern[i] = 1;
			}
		}
		pattern = as.integer(pattern);
		result = list(pattern=pattern, hypotheses=hypotheses, function.name=function.name, function.inputs=function.inputs, fisher.pvalues=fisher.pvalues);
		return(result);
	}
	else {
		stop("[ERR] Either the genotypes or the pattern not provided! No hypothesis will be created.");
	}
	return(NA)
}

# OR hypothesis
#' @export
"OR" <-
function( ... ) {
	# look for the global variables named lifting.genotypes and lifting.annotations
	genotypes = lifting.genotypes;
	annotations = lifting.annotations;
	if(!is.null(genotypes) && !is.null(annotations) && length(list(...))>0) {
		# get the vector of the clauses of the pattern from the genotypes
		result = aux.log(genotypes,annotations,"OR",...);
		curr_genotypes = result$curr_genotypes;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		# save the fisher exact tests
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				curr.p.value = 1 - curr.test[1];
	
				fisher.pvalues = append(fisher.pvalues,curr.p.value)
			}
			assign("fisher.pvalues", fisher.pvalues, envir=.GlobalEnv)
		}
		
		# evaluate the OR operator
		pattern = rep(0,nrow(genotypes));
		for (i in 1:nrow(genotypes)) {
			pattern[i] = sum(curr_genotypes[i,]);
			if(pattern[i]>1) {
				pattern[i] = 1;
			}
		}
		pattern = as.integer(pattern);
		result = list(pattern=pattern, hypotheses=hypotheses, function.name=function.name, function.inputs=function.inputs, fisher.pvalues=fisher.pvalues);
		return(result);
	}
	else {
		stop("[ERR] Either the genotypes or the pattern not provided! No hypothesis will be created.");
	}
	return(NA);
}

# XOR hypothesis
#' @export
"XOR" <-
function( ... ) {
	#look for the global variables named lifting.genotypes and lifting.annotations
	genotypes = lifting.genotypes;
	annotations = lifting.annotations;
	if(!is.null(genotypes) && !is.null(annotations) && length(list(...))>0) {
		#get the vector of the clauses of the pattern from the genotypes
		result = aux.log(genotypes,annotations,"XOR",...);
		curr_genotypes = result$curr_genotypes;
		hypotheses = result$hypotheses;
		function.name = result$function.name;
		function.inputs = result$function.inputs;
		
		# save the fisher exact tests
		if(length(result$tests)>0) {
			for(i in 1:length(result$tests)) {
				curr.test = result$tests[[i]];
				odds.ratio = curr.test[2]
		
				if(curr.test[2]>=0) {
					curr.pvalue = 1;
				}
				else {
					curr.pvalue = curr.test[1];
				}
				
				fisher.pvalues = append(fisher.pvalues,curr.pvalue)
			}
			assign("fisher.pvalues", fisher.pvalues, envir=.GlobalEnv)
		}
		
		# evaluate the XOR operator
		pattern = rep(0,nrow(genotypes));
		for (i in 1:nrow(genotypes)) {
			pattern[i] = sum(curr_genotypes[i,]);
			if(pattern[i]>1) {
				pattern[i] = 0;
			}
		}
		pattern = as.integer(pattern);
		result = list(pattern=pattern, hypotheses=hypotheses, function.name=function.name, function.inputs=function.inputs, fisher.pvalues=fisher.pvalues);
		return(result);
	}
	else {
		stop("[ERR] Either the genotypes or the pattern not provided! No hypothesis will be created.");
	}
	return(NA);
}

#### end of file -- lifting.operators.R

#### emap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the position in the genotypes of the event referring to the given label.event
emap = function( label.event, genotypes, annotations ) {
	col.num = -1;
	events.name = "";
	if(!is.null(genotypes) && !is.null(annotations)) {
		if(label.event[2]!="*") {
			curr.events = which(annotations[,"event"]==label.event[1]&annotations[,"type"]==label.event[2]);
		}
		else {
			curr.events = which(annotations[,"event"]==label.event[1]);
		}
		if(length(curr.events)>0) {
			events.name = names(curr.events);
			col.num = which(colnames(genotypes)%in%events.name);
		}
	}
	else {
		stop("[ERR] A genotypes must be available in order to define any hypothesis!",call.=FALSE);
	}
	results = list(col.num=col.num,events.name=events.name);
	return(results);
}

#### end of file -- emap.R

#### get.lifted.pattern.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the adjacency matrix of the pattern given the list of edges involving it
"get.lifted.pattern" <-
function( lifted.edges ) {
	#structure to save the adjacency matrix
	lifted.adj.matrix = array(0,c(length(unique(c(lifted.edges[,1],lifted.edges[,2]))),length(unique(c(lifted.edges[,1],lifted.edges[,2])))));
	rownames(lifted.adj.matrix) = unique(c(lifted.edges[,1],lifted.edges[,2]));
	colnames(lifted.adj.matrix) = rownames(lifted.adj.matrix);
	#build the matrix given the lifted.edges
	for(i in 1:nrow(lifted.edges)) {
		lifted.adj.matrix[lifted.edges[i,1],lifted.edges[i,2]] = 1;
	}
	return(lifted.adj.matrix);
}

#### end of file -- get.lifted.pattern.R

#### hypothesis.cycles.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# evaluate cycles involving any hypothesis
# INPUT:
# data: input genotypes and its hypotheses
# adj.matrix: adjacency matrix of the reconstructed topology
# hypotheses.labels: label of all the existing hypotheses
# weights.matrix: weights of any edge in the topology
"hypothesis.evaluate.cycles" <-
function( data, adj.matrix, hypotheses.labels, weights.matrix ) {
	
	# create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
	# create a map structure where to save the atomic events of each hypothesis
	matomic = new.env(hash=TRUE,parent=emptyenv());
	
	# create a map structure where to save the hypotheses of each atomic event
	mhypotheses = new.env(hash=TRUE,parent=emptyenv());
	
	# evaluate all the existing hypotheses
	for (i in 1:length(hypotheses.labels)) {
		
		#evaluate the current hypothesis
		connections = hypothesis.connections(adj.matrix, hypotheses.labels[i]);
		connections = hypothesis.expand.connections(label=hypotheses.labels[i],events=pattern.events(data,hypotheses.labels[i]),incoming=connections$incoming,outgoing=connections$outgoing,hnames=colnames(adj.matrix),matomic=matomic,weights.matrix=weights.matrix);
		
		#save the results for the current hypothesis
		ordered.weights = c(ordered.weights,connections$ordered.weights);
		ordered.edges = c(ordered.edges,connections$ordered.edges);
		matomic = connections$matomic;
		
	}
	
	#add to the map the link between atomic to pattern
	for (i in 1:ncol(adj.matrix)) {
		if(!is.null(data$atoms[[colnames(adj.matrix)[i]]])) {
			#add to the map the hypotheses of this atomic event
			mhypotheses[[toString(i)]] = which(colnames(adj.matrix)%in%data$atoms[[colnames(adj.matrix)[i]]]);
		}
	}
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,matomic=matomic,mhypotheses=mhypotheses));
}

# given the adj.matrix, return the incoming and outgoing connections for any hypothesis
# INPUT:
# adj.matrix: adjacency matrix of the topology
# hypotheses.label: label of the hypothesis
"hypothesis.connections" <-
function( adj.matrix, hypotheses.label ) {
	
	hypotheses.label = hypotheses.label[hypotheses.label %in% rownames(adj.matrix)]
	
	incoming = rownames(adj.matrix)[which(adj.matrix[,hypotheses.label]==1)];
	outgoing = colnames(adj.matrix)[which(adj.matrix[hypotheses.label,]==1)];
	connections = list(incoming=incoming,outgoing=outgoing);
	
	return(connections);

}

# expand and enumerate all the connections incoming or outgoing an hypothesis
# INPUT:
# label: name of the hypothesis
# events: events in the hypothesis
# incoming: incoming connections
# outgoing: outgoing connections
# weights.matrix: weights of any edge in the topology
"hypothesis.expand.connections" <-
function(label, events, incoming, outgoing, hnames, matomic, weights.matrix) {
	
	# create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
    # get the position of the hypothesis
    hypothesis.pos = which(hnames==label);

	# evalutate the incoming and outgoing connections
    curr.edge.pos = 0;
    if(length(incoming)>0) {
		for(i in 1:length(incoming)) {
			ordered.weights = rbind(ordered.weights,weights.matrix[which(hnames==incoming[i]),hypothesis.pos]);
			curr.edge.pos = curr.edge.pos + 1;
			new.edge <- array(0, c(2,1));
        		new.edge[1,1] = which(hnames==incoming[i]);
        		new.edge[2,1] = hypothesis.pos;
        		ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
    if(length(outgoing)>0) {
		for(i in 1:length(outgoing)) {
			ordered.weights = rbind(ordered.weights,weights.matrix[hypothesis.pos,which(hnames==outgoing[i])]);
			curr.edge.pos = curr.edge.pos + 1;
        		new.edge <- array(0, c(2,1));
        		new.edge[1,1] = hypothesis.pos;
        		new.edge[2,1] = which(hnames==outgoing[i]);
        		ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
	
	if(length(hypothesis.pos) > 0)
		matomic[[toString(hypothesis.pos)]] = which(hnames%in%events)	#add to the map the atomic events of this hypothesis
	else
		print('hypothesis.pos == 0!')
	
	# return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,matomic=matomic));
	
}

# given the hypotheses and the adj.matrix, return the updated adj.matrix
"hypothesis.adj.matrix" <-
function(hypotheses, adj.matrix) {
	
	if(!is.na(hypotheses[1])) {
		
		# set the invalid entries in the adj.matrix
		# hypotheses can not be causing other hypotheses
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix),(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix)] = 0;
		
		# consider the given hypotheses only against the specified possible effects
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix),1:(ncol(adj.matrix)-hypotheses$num.hypotheses)] = 0
		adj.matrix[1:(ncol(adj.matrix)-hypotheses$num.hypotheses),(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix)] = 0
		
		# set the elements from the hlist
		for (i in 1:nrow(hypotheses$hlist)) {
			cause = which(colnames(adj.matrix)%in%hypotheses$hlist[i,"cause"]);
			effect = which(colnames(adj.matrix)%in%hypotheses$hlist[i,"effect"]);
			if(length(cause)>0 && length(effect)>0) {
				adj.matrix[cause,effect] = 1;
			}
		}
	}
	return(adj.matrix);
	
}


#### end of file -- hypothesis.cycles.R

testing = function(data, g1, g2) {

	# Dataframe di tutto il genotypes
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
	
	cat('genotypes\n')
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
					
				}
			}
		}
	}
	
	return(results)
	
}
