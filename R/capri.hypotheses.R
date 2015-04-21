#### emap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the position of the event of a given label.event
emap = function( label.event, dataset, annotations ) {
	col.num = -1;
	events.name = "";
	if(!is.null(dataset) && !is.null(annotations)) {
		if(label.event[2]!="*") {
			curr.events = which(annotations[,"event"]==label.event[1] & annotations[,"type"]==label.event[2]);
		}
		else {
			curr.events = which(annotations[,"event"]==label.event[1]);
		}
		if(length(curr.events)>0) {
			events.name = names(curr.events);
			col.num = which(colnames(dataset)%in%events.name);
		}
	}
	else {
		stop("[ERR] The dataset must be available to define hypotheses!",call.=FALSE);
	}
	results = list(col.num=col.num,events.name=events.name);
	return(results);
}

#### end of file -- emap.R

#### get.lifted.formula.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the adjacency matrix of the formula given the list of edges
"get.lifted.formula" <-
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

#### end of file -- get.lifted.formula.R


#### hypothesis.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Add a new hypothesis by creating a new causal event and adding it to the dateset
#' @export
"hypothesis.add" <-
  # function( data, label.formula, pvalue = 0.05, lifted.formula,  ... ) {
  function( data, label.formula, lifted.formula,  ... ) {

	pvalue = 0.05
	
    label.effect = list(...);
    #print(label.effect)
    
    # save the needed data structures
    if(!is.null(data$genotypes) && !is.null(data$annotations)) {
      dataset = data$genotypes;
      annotations = data$annotations;
    }
    else {
      dataset = NULL;
      annotations = NULL;
    }
    if(!is.null(data$hypotheses)) {
      hypotheses = data$hypotheses;
    }
    else {
      hypotheses = NA;
    }
    
    # add the hypothesis only if all the inputs are correctly provided
    if(!is.null(dataset) && !is.null(annotations)) {
      # the Boolean functions look for a global variable named lifting.dataset
      # if there are already global variables named as the ones used here, make the backup of them
      do.roll.back.lifting.dataset = FALSE;
      do.roll.back.lifting.annotations = FALSE;
      do.roll.back.lifting.edges = FALSE;
      
      # I need a global variable to save the dataset of the lifted formula
      # if there is already a global variable named lifting.dataset, make the backup of it
      if(exists("lifting.dataset")) {
        roll.back.lifting.dataset = lifting.dataset;
        do.roll.back.lifting.dataset = TRUE;
      }
      assign("lifting.dataset",dataset,envir=.GlobalEnv);
      
      # I need a global variable to save the annotations of the lifted formula
      # if there is already a global variable named lifting.annotations, make the backup of it
      if(exists("lifting.annotations")) {
        roll.back.lifting.annotations = lifting.annotations;
        do.roll.back.lifting.annotations = TRUE;
      }
      assign("lifting.annotations",annotations,envir=.GlobalEnv);
      
      # I need a global variable to save the edges of the lifted formula
      # if there is already a global variable named lifting.edges, make the backup of it
      
      do.roll.back.lifting.dataset = FALSE; # <- ???? why????
      
      if(exists("lifting.edges")) {
        roll.back.lifting.edges = lifting.edges;
        do.roll.back.lifting.edges = TRUE;
      }
      assign("lifting.edges",NULL,envir=.GlobalEnv);
      
      # I need a global variable to save the pvalue to be used for the fisher exact test
      # if there is already a global variable named lifting.pvalue, make the backup of it
      do.roll.back.lifting.pvalue = FALSE;      
      if(exists("lifting.pvalue")) {
        roll.back.lifting.pvalue = lifting.pvalue;
        do.roll.back.lifting.pvalue = TRUE;
      }
      assign("lifting.pvalue",pvalue,envir=.GlobalEnv);
      
      #print('lifted.formula')
      #print(lifted.formula)
      
      ## test
      #lifted.formula = eval(lifted.formula)
      
      
      # save the lifted dataset and its hypotheses for the current formula
      curr_formula = lifted.formula$formula;
      curr_hypotheses = lifted.formula$hypotheses;


      #print('cur formula')
      #print(curr_formula)
      #print('cur_hypo')
      #print(curr_hypotheses)
      
      # save the edges of the lifted formula
      hstructure = lifting.edges;
      
      # roll back to the previous value of the global variable lifting.dataset if any or remove it
      if(do.roll.back.lifting.dataset) {
        assign("lifting.dataset",roll.back.lifting.dataset,envir=.GlobalEnv);
      }
      else {
        rm(lifting.dataset,pos=".GlobalEnv");
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
      # roll back to the previous value of the global variable lifting.pvalue if any or remove it
      if(do.roll.back.lifting.pvalue) {
        assign("lifting.pvalue",roll.back.lifting.pvalue,envir=.GlobalEnv);
      }
      else {
        rm(lifting.pvalue,pos=".GlobalEnv");
      }
      
      # set the hypotheses number
      if(!is.na(hypotheses[1])) {
        num.hypotheses = hypotheses$num.hypotheses;
      }
      else {
        num.hypotheses = 0;
      }
      # * is a special label.effect which indicates to use all the events as effects for this formula
      is.to.all.effects = FALSE;
            
      if(label.effect[[1]][1]=="*") {
        label.effect = colnames(dataset)[1:(length(colnames(dataset))-num.hypotheses)];
        #any event can not be both causes and effects for the formula to be well-formed
        label.effect = list(label.effect[-which((label.effect%in%unlist(curr_hypotheses$llist)))]);
        is.to.all.effects = TRUE;
        
        if(length(label.effect)==0) 
        {
          stop(paste("[ERR] Missing list of effects to test or wildcard \'*\'.", sep=''));
        }
      }
      # check the formula to be well-formed
      all.col.nums = vector();
      if(length(label.effect)==0) {
        stop(paste("[ERR] Missing list of effects or wildcard \'*\'.", sep=''));
      }
      else {
        #check the effects of the formula to be well-formed
        for (i in 1:length(label.effect)) {
          curr.label.effect = label.effect[[i]];
          if(is.to.all.effects==FALSE) {
            col.num = -1;
            if(length(curr.label.effect)==1) {
              event.map = emap(c(curr.label.effect,"*"),dataset,annotations);
              col.num = event.map$col.num;
              events.name = event.map$events.name;
            }
            else if(length(curr.label.effect)==2) {
              event.map = emap(curr.label.effect,dataset,annotations);
              col.num = event.map$col.num;
              events.name = event.map$events.name;
            }
          }
          else {
            col.num = which(colnames(dataset)%in%curr.label.effect);
            if(length(col.num)==0) {
              col.num = -1;
            }
            events.name = curr.label.effect;
          }
          #check the effect to be a valid event
          if(col.num[1]==-1) {            
            stop(paste("[ERR] Unknown gene among effects: \"", curr.label.effect,
                       "\".",sep=''));
          }
          all.col.nums = append(all.col.nums,col.num);
          #check the formula to be well-formed
          #if the effect is in the formula, the formula is not well-formed
          if(length(which(unlist(curr_hypotheses$llist)%in%events.name))>0) {
                  stop(paste("[ERR] Bad formed formula, event \"", curr.label.effect,
                       "\" yields a loop.",,sep=''));          
            }
        }
      }
      
      # look for duplicated effects in the formula
      if(anyDuplicated(all.col.nums)>0) 
        {
        stop(paste("[ERR] Bad formed formula, duplicated events ", 
                   paste(label.effect[duplicated(label.effect)], collapse=', ', sep=''),
                   "within effects.", sep=''));          
        }
      #check that the we are not duplicating any name by adding the new hypothesis
      if(length(which(colnames(dataset)==label.formula))>0) 
      {
        stop(paste("[ERR] This hypothesis already exists.", sep=''));
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
      if(marginal.probs[ncol(dataset)]==0 || marginal.probs[ncol(dataset)]==1) 
        {        
        stop(paste("[ERR] Formula has marginal probability ", marginal.probs[ncol(dataset)], 
                   ", but should be in (0,1).", sep=''));
      }
      #check that the formula does not duplicate any existing column
      i = ncol(dataset);
      for(j in 1:ncol(dataset)) {
        #if the edge is valid, i.e., not self cause
        if(i!=j) {
          #if the two considered events are not distinguishable
          if((joint.probs[i,j]/marginal.probs[i])==1 && (joint.probs[i,j]/marginal.probs[j])==1) 
          {
            stop(paste("[ERR] Pattern duplicates ", paste(as.events(data)[j, ], collapse=' ', sep=''), 
                       ".", sep=''));
          }
        }
      }
      #now I can finally add the hypothesis
      colnames(dataset)[ncol(dataset)] = label.formula;
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
        curr.label.effect = label.effect[[i]];
        if(is.to.all.effects==FALSE) {
          if(length(curr.label.effect)==1) {
            event.map = emap(c(curr.label.effect,"*"),dataset,annotations);
            col.num = event.map$col.num;
          }
          else if(length(curr.label.effect)==2) {
            event.map = emap(curr.label.effect,dataset,annotations);
            col.num = event.map$col.num;
          }
        }
        else {
          col.num = which(colnames(dataset)%in%curr.label.effect);
          if(length(col.num)==0) {
            col.num = -1;
          }
        }
        for (j in 1:length(col.num)) {
          hypotheses$hlist = rbind(hypotheses$hlist,t(c(colnames(dataset)[ncol(dataset)],colnames(dataset)[col.num[j]])));
        }
        if(is.null(colnames(hypotheses$hlist))) {
          colnames(hypotheses$hlist) = c("cause","effect");
        }
      }
      #create the list of hypotheses' structures
      if(length(hypotheses$hstructure)==0) {
        hypotheses$hstructure = new.env(hash=TRUE,parent=emptyenv());
      }
      
      #add the atoms in the hypothesis
      if(length(hypotheses$patterns)==0) {
      	hypotheses$patterns = list()
      }
      hypotheses$patterns[label.formula] = lifted.formula$hypotheses$llist;
      
      hypotheses$hstructure[[label.formula]] = get.lifted.formula(hstructure);
      #add the new hypothesis in the annotations
      annotations = rbind(data$annotations,c("Hypothesis", label.formula));
      rownames(annotations)[nrow(annotations)] = label.formula;
      #add the color of the type "Hypothesis" is not already defined
      if(any(rownames(data$types)=="Hypothesis")==FALSE) {
        types = rbind(data$types, 'slateblue');
        rownames(types)[nrow(types)] = "Hypothesis";
        data$types = types;
      }
      
      #add the hypotheses in the atoms
      if(length(hypotheses$atoms)==0) {
      	hypotheses$atoms = vector(mode="list",length=(ncol(dataset)-hypotheses$num.hypotheses));
      	names(hypotheses$atoms) = colnames(dataset)[1:(ncol(dataset)-hypotheses$num.hypotheses)];
      }
      atoms.in.formula = which(names(hypotheses$atoms)%in%unlist(hypotheses$patterns[label.formula]));
      if(length(atoms.in.formula)>0) {
      	for (i in 1:length(atoms.in.formula)) {
      		hypotheses$atoms[[atoms.in.formula[i]]] = append(hypotheses$atoms[[atoms.in.formula[i]]], label.formula);
      	}
      }
      
      
      #return the new data as result
      data$genotypes = dataset;
      data$hypotheses = hypotheses;
      data$annotations = annotations;
      return(data);
    }
    else {
      stop("[ERR] Missing dataset or formula.");
    }
    return(NA);
  }

#### end of file -- hypothesis.add.R

# commenti

#' @export
hypothesis.add.group = function(x, FUN, group, dim.min = 2, dim.max = length(group), min.prob = 0, ...) {
	op = deparse(substitute(FUN))

	#print(length(unlist(group)))
	#print(unlist(group))

	effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste(effect, collapse = ", ")

	
	cat("*** Adding Group Hypotheses\n")
	cat('Group:', paste(group, collapse = ", ", sep = ""))
	cat(' Function:', op)
	cat(' Effect:', effect, '\n')
	flush.console()
	
	# group %in% as.events(x)[, 'event']


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
		warning("No hypothesis can be created for groups with less than 2 elements.")
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

	# Get an analytical formula... !
	tot.patterns = 0
	for (i in min.groupsize:max.groupsize) tot.patterns = tot.patterns + ncol(combn(unlist(group), i))
	
	# create a progress bar
	cat('Generating ', tot.patterns ,'patterns [size: min =', max.groupsize,' -  max =', max.groupsize, '].\n')
		
	# pb <- txtProgressBar(0, tot.patterns, style = 3)
	flush.console()

	pbPos = 0
	for (i in min.groupsize:max.groupsize) {
		gr = combn(unlist(group), i)
	
		# print(gr)
		
		
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

				hypo.add = paste0("hypothesis.add(x, label.formula = '", op, "_", hypo.name, "', lifted.formula = ", op, "(", hypo.genes, "), ", effect, ")")


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
hypothesis.add.homologous = function(x, ..., genes = as.genes(x), FUN = "OR") {
	# in questa funzione, per ogni gene che ha più di un tipo di alterazione
	# aggiungo l'OR

	effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste(effect, collapse = ", ")

	hom.group = lapply(genes, function(g, x) {
		if (nevents(x, genes = g) > 1) 
			T
		else F
	}, x)
	hom.group = genes[unlist(hom.group)]

	cat("*** Adding hyoptheses for Homolgous Patterns\n")
	cat('Genes:', paste(hom.group, collapse = ", ", sep = ""))
	cat(' Function:', FUN)
	cat(' Effect:', effect, '\n')
	flush.console()
	
	# create a progress bar
	pb <- txtProgressBar(0, length(hom.group), style = 3)

	error.summary = data.frame()

	for (i in 1:length(hom.group)) {

		#start the progress bar
		setTxtProgressBar(pb, i)

		# Check if the joint probability of homologous events is > 0, if
		# yes the event will be added as 'OR', otherwise 'XOR'
		if( any(
				rowSums(as.gene(x, genes = hom.group[[i]])) > 1) 
			)
		FUN = 'OR'
		else FUN = 'XOR'				

		hypo.add = paste0("hypothesis.add(x, label.formula = '", FUN, "_", hom.group[[i]], "', lifted.formula = ", FUN, "('", hom.group[[i]], "'), ", effect, 
			")")

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


#### hypothesis.cycles.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#evaluate cycles involving any hypothesis
#INPUT:
#data: input dataset and its hypotheses
#adj.matrix: adjacency matrix of the topology
#hypotheses.labels: label of the existing hypotheses
"hypothesis.evaluate.cycles" <-
function(data, adj.matrix, hypotheses.labels, weights.matrix) {
	
	#create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
	#create a map structure where to save the atomic events of each hypothesis
	matomic = new.env(hash=TRUE,parent=emptyenv());
	
	#create a map structure where to save the hypotheses of each atomic event
	mhypotheses = new.env(hash=TRUE,parent=emptyenv());
	
	#evaluate all the existing hypotheses
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

#given the adj.matrix, return the incoming and outgoing connections for any hypothesis
#INPUT:
#adj.matrix: adjacency matrix of the topology
#hypotheses.label: label of the hypothesis
"hypothesis.connections" <-
function(adj.matrix, hypotheses.label) {
	# cat('\nhl', hypotheses.label, '\n nomi colonne:\n',
	 # paste(colnames(adj.matrix)), 'asd')
	
	#### EROS FIX
	# print('*** PRE')
	# print(hypotheses.label)
	# print('*** POST')
	# foo = hypotheses.label
	hypotheses.label = hypotheses.label[hypotheses.label %in% rownames(adj.matrix)]
	# print(hypotheses.label)
	# if(length(hypotheses.label) == 0) {
		# print(rownames(adj.matrix))
		# print(foo)
		# print(foo %in% rownames(adj.matrix))
		# }
	
	
	incoming = rownames(adj.matrix)[which(adj.matrix[,hypotheses.label]==1)];
	outgoing = colnames(adj.matrix)[which(adj.matrix[hypotheses.label,]==1)];
	connections = list(incoming=incoming,outgoing=outgoing);
	
	#### EROS FIX
	# print(connections)
	return(connections);
}

#expand and enumerate all the connections incoming or outgoing an hypothesis
#INPUT:
#label: name of the hypothesis
#events: events in the hypothesis
#incoming: incoming connections
#outgoing: outgoing connections
"hypothesis.expand.connections" <-
function(label, events, incoming, outgoing, hnames, matomic, weights.matrix) {
	
	#create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
    #get the position of the hypothesis
    hypothesis.pos = which(hnames==label);

	#evalutate the incoming and outgoing connections
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
	
	#### EROS FIX
	# print('****')
	# print(hypothesis.pos)
	# print('hnames')
	# print(hnames)
	# print('events')
	# print(events)
	if(length(hypothesis.pos) > 0)
		matomic[[toString(hypothesis.pos)]] = which(hnames%in%events)	#add to the map the atomic events of this hypothesis
	else
		print('hypothesis.pos == 0!')
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,matomic=matomic));
}

#given the hypotheses and the adj.matrix, return the updated adj.matrix
"hypothesis.adj.matrix" <-
function(hypotheses, adj.matrix) {
	if(!is.na(hypotheses[1])) {
		# set the invalid entries in the adj.matrix
		# hypotheses can not be causing other hypotheses
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix),(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix)] = 0;
		# consider the given hypotheses only toward the specified possible effects
		hypotheses.matrix = array(0,c(hypotheses$num.hypotheses,ncol(adj.matrix)-hypotheses$num.hypotheses));		
		for (i in 1:nrow(hypotheses$hlist)) {
			cause = which(hypotheses$hlist[i,1]==colnames(adj.matrix));
			effect = which(hypotheses$hlist[i,2]==colnames(adj.matrix));
			if(length(cause)>0 && length(effect)>0) {
				hypotheses.matrix[cause-ncol(adj.matrix)+hypotheses$num.hypotheses,effect] = 1;
			}
		}
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix),1:(ncol(adj.matrix)-hypotheses$num.hypotheses)] = hypotheses.matrix;
		for(j in (ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix)) {
			for(k in 1:(ncol(adj.matrix)-hypotheses$num.hypotheses)) {
				if(adj.matrix[j,k] == 0) {
					adj.matrix[k,j] = 0;
				}
			}
		}
	}
	return(adj.matrix);
}


#### end of file -- hypothesis.cycles.R


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
#' @export
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
#' @export
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
#' @export
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

#' @import igraph
hypotheses.expansion <- function(input_matrix, 
                                 map = list(),
                                 hidden_and = T,
                                 expand = T,
                                 events = NULL,
                                 conf_matrix = NULL,
                                 skip.disconnected = TRUE
                                 ) {
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  num_hypos = length(map)
  
  # get node list
  node_list <- colnames(input_matrix)
  #print('input matrix')
  #print(input_matrix)
  
  # cut input matrix
  margin = length(node_list) - num_hypos
  hypos_new_name = list()
  
  # da finire!!!
  if(is.vector(events) && F) {
    cat('\n remove event is broken!!!\n')
    print(events)
    min_graph = graph.adjacency(input_matrix)
    graph <- igraph.to.graphNEL(min_graph)
    edge_names = edgeNames(graph)
    print(edge_names)
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      print(edge)
      from = edge[1]
      to = edge[2]
      check_from = any(unlist(strsplit(from, '_')) %in% events)
      check_to = any(unlist(strsplit(to, '_')) %in% events)
      if (!(check_from && check_to)) {
        input_matrix[from, to] = 0
      }
    }
    print(input_matrix)
  }
  
  

  cat('*** Hypos expansion:')
  # check if there are hypotheses
  if (num_hypos == 0 || !expand) {
    # if no hypos do nothings..
    min_matrix = input_matrix
  } else {
    # ..else expand them
    min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]
    
    # create graph from matrix
    min_graph = graph.adjacency(min_matrix)
    
    
    # foreach hypothesis
    # print(ls(map))
    # for (h in ls(map)) {
    # print(h)
    # }
    
    
    for (h in ls(map)) {
      

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

      # change names in confidence matrix according to hypotesis
      if(!is.null(conf_matrix)) {
        rownames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
        colnames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
      }

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

        # edge to reconstruct
        h_edge <- input_matrix[, h]
        initial_node_up <- names(h_edge)[which(h_edge==1)]

        # add this graph to main graph
        min_graph = graph.union(min_graph, hypo_graph_pre)

        # recreate lost edge
        for (node in initial_node_up) {
          min_graph <- min_graph + edge(node, final_node)
        }

      }

      
      # display down hypo
     # if (display.up || display.down) {
     if (display.down) {
     

        # edge to reconstruct
        h_edge <- input_matrix[h,]
        final_node <- names(h_edge)[which(h_edge==1)]

        # add this graph to main graph
        min_graph = graph.union(min_graph, hypo_graph)
        

      }

      # reconnect down hypo
      if (display.down) {
        # print(final_node)
        # recreate lost edge
        for (node in final_node) {
          min_graph <- min_graph + edge(initial_node, node)
        }
      }
      
    }
    min_matrix = get.adjacency(min_graph, sparse = F)
    
  }
  
  cat(' done')

  
  
  # now expand the hidden AND
  #print(min_matrix)
  
  if(hidden_and == F) {
    # sort col and row (igraph wants the same order)
    min_matrix = min_matrix[,order(colnames(min_matrix))]
    min_matrix = min_matrix[order(rownames(min_matrix)),]
    
    # print(min_matrix)
    if(!is.null(conf_matrix)) {
      return(list(min_matrix, hypos_new_name, conf_matrix))
    }
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
  
  cat(' done')
  
  # sort col and row (igraph wants the same order)
  and_matrix = and_matrix[,order(colnames(and_matrix))]
  and_matrix = and_matrix[order(rownames(and_matrix)),]
  
  # print(and_matrix)
  if(!is.null(conf_matrix)) {
    return(list(and_matrix, hypos_new_name, conf_matrix))
  }
  return(list(and_matrix, hypos_new_name))
}

