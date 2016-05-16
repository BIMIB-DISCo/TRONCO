#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# reconstruct the best topology based on probabilistic causation and Prim algorithm
# @title prim.fit
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent should I be verbose?
# @return topology: the reconstructed tree topology
#
prim.fit <- function(dataset,
                      regularization = "no_reg",
                      do.boot = TRUE,
                      nboot = 100,
                      pvalue = 0.05,
                      min.boot = 3,
                      min.stat = TRUE,
                      boot.seed = NULL,
                      silent = FALSE ) {

    ## Start the clock to measure the execution time.
    
    ptm = proc.time();

    ## Structure with the set of valid edges
    ## I start from the complete graph, i.e., I have no prior and all
    ## the connections are possibly causal.
    
    adj.matrix = array(1, c(ncol(dataset), ncol(dataset)));
    colnames(adj.matrix) = colnames(dataset);
    rownames(adj.matrix) = colnames(dataset);

    ## The diagonal of the adjacency matrix should not be considered,
    ## i.e., no self cause is allowed.
    
    diag(adj.matrix) = 0;

    ## Check if the dataset is valid.
    
    valid.dataset = check.dataset(dataset, adj.matrix, FALSE);
    adj.matrix = valid.dataset$adj.matrix;
    invalid.events = valid.dataset$invalid.events;

    ## Reconstruct the prima facie topology
    ## Should I perform bootstrap? Yes if TRUE, no otherwise.
    
    if (do.boot == TRUE) {
        if (!silent)
            cat('*** Bootstraping selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.do.boot(dataset,
                                            NA,
                                            nboot,
                                            pvalue,
                                            adj.matrix,
                                            min.boot,
                                            min.stat,
                                            boot.seed,
                                            silent);
    } else {
        if (!silent)
            cat('*** Computing selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.no.boot(dataset,
                                            NA,
                                            adj.matrix,
                                            silent);
    }

    ## Add back in any connection invalid for the probability raising
    ## theory.
    
    if (length(invalid.events) > 0) {
        for (i in 1:nrow(invalid.events)) {
            prima.facie.parents$adj.matrix$adj.matrix.acyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1;
            prima.facie.parents$adj.matrix$adj.matrix.cyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1;
        }
    }
    adj.matrix.prima.facie =
        prima.facie.parents$adj.matrix$adj.matrix.acyclic

    ## Perform the likelihood fit with the required strategy.
    
    model = list();
    for (reg in regularization) {

        ## Perform the likelihood fit with the chosen regularization
        ## score on the prima facie topology.
        
        if (!silent) {
            cat('*** Performing likelihood-fit with regularization:',reg,'.\n')
        }
        best.parents =
            perform.likelihood.fit.prim(dataset,
                adj.matrix.prima.facie,
                regularization = reg,
                marginal.probs = prima.facie.parents$marginal.probs,
                joint.probs = prima.facie.parents$joint.probs)

        ## Set the structure to save the conditional probabilities of
        ## the reconstructed topology.

        reconstructed.model = create.model(dataset,
            best.parents,
            prima.facie.parents)

        model.name = paste('prim', reg, sep='_')
        model[[model.name]] = reconstructed.model
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "PRIM",
             regularization = regularization,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent)

    ## Return the results.
    
    topology =
        list(dataset = dataset,
             hypotheses = NA,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             parameters = parameters,
             execution.time = (proc.time() - ptm))
    topology = rename.reconstruction.fields(topology, dataset)
    return(topology)
}


# reconstruct the best causal topology by Prim algorithm combined with probabilistic causation
# title perform.likelihood.fit.prim
# param dataset a valid dataset
# param adj.matrix the adjacency matrix of the prima facie causes
# param regularization regularization term to be used in the likelihood fit
# param command type of search, either hill climbing (hc) or tabu (tabu)
# return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.prim = function(dataset,
                                       adj.matrix,
                                       regularization,
                                       command = "hc",
                                       marginal.probs,
                                       joint.probs){

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
       
    # create the graph of valid edges
    curr.valid.adj.matrix = array(0, c(nrow(adj.matrix), ncol(adj.matrix)))
    cont = 0
    for (i in 1:nrow(curr.valid.adj.matrix)) {
        for (j in i:nrow(curr.valid.adj.matrix)) {
        	# we want an undirected prima facie graph, 
        	# so we consider all the edges where there is a directed arc
            if (!(adj.matrix[i,j] == 0 && adj.matrix[j,i] == 0)) {
                cont = cont + 1
                curr.valid.adj.matrix[i,j] = 1
                curr.valid.adj.matrix[j,i] = 1
            }
        }
    }
    
    # if the graph has at least one edge, build to weighted igraph object
    if (cont > 0) {
        curr.graph = graph.adjacency(curr.valid.adj.matrix, mode = "undirected")
        all_edges = get.edgelist(curr.graph)
        # set the weights to the edges
        new_weights = NULL
        
        if (cont == 1) {
        	# consider the current arc
        	i = all_edges[1,1]
        	j = all_edges[1,2]
	        # if the event is valid
	        if(joint.probs[i,j]>=0) {
	            new_score = compute.mi.score(joint.probs[i,j],marginal.probs[i],marginal.probs[j]) # log(joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j]))
	        }
	        # else, if the two events are indistinguishable
	        # put the higher score
	        else if(joint.probs[i,j]<0) {
	            new_score = 1 # Inf
	        }
            new_weights = new_score # mutinformation(data[ ,all_edges[1,1]], data[ ,all_edges[1,2]])
        } else {
            for (i in 1:nrow(all_edges)) {
            	# consider the current arc
        	    curr_i = all_edges[i,1]
        	    curr_j = all_edges[i,2]
		        # if the event is valid
		        if(joint.probs[curr_i,curr_j]>=0) {
		            new_score = compute.mi.score(joint.probs[curr_i,curr_j],marginal.probs[curr_i],marginal.probs[curr_j]) # log(joint.probs[curr_i,curr_j]/(marginal.probs[curr_i]*marginal.probs[curr_j]))
		        }
		        # else, if the two events are indistinguishable
		        # put the higher score
		        else if(joint.probs[curr_i,curr_j]<0) {
		            new_score = 1 # Inf
		        }
                new_weights = c(new_weights,new_score) # mutinformation(data[ ,all_edges[i,1]], data[ ,all_edges[i,2]]))
            }
        }
        
        # # set the weights to the graph
        # if(length(new_weights[new_weights!=Inf])>0) {
            # inf.scores = which(new_weights==Inf)
            # max_score = max(new_weights[new_weights!=Inf])
            # prim_scores = (max_score - new_weights) / max_score
            # prim_scores[inf.scores] = 0
            # E(curr.graph)$weight = prim_scores
        # }
        # else {
        	# inf.scores = which(new_weights==Inf)
        	# new_weights[inf.scores] = 0
            # E(curr.graph)$weight = new_weights
        # }
        
        # set the weights to the graph
        E(curr.graph)$weight = 1 - new_weights # max(new_weights) - new_weights
        
        # get the minimum spanning tree by Prim algorithm
        curr.valid.adj.matrix = as.matrix(get.adjacency(minimum.spanning.tree(curr.graph, 
                                                                              algorithm="prim")))

        # build the matrix of the priors
        new_prior_matrix = array(0, c(nrow(adj.matrix), ncol(adj.matrix)))
        rownames(new_prior_matrix) = colnames(dataset)
        colnames(new_prior_matrix) = colnames(dataset)
        for (i in 1:nrow(curr.valid.adj.matrix)) {
            for (j in i:ncol(curr.valid.adj.matrix)) {
                if (i != j) {
                    if (curr.valid.adj.matrix[i,j] == 1) {
                        if (adj.matrix[i,j] == 1) {
                            new_prior_matrix[i,j] = 1
                        } else if(adj.matrix[j,i] == 1) {
                            new_prior_matrix[j,i] = 1
                        }
                    }
                }
            }
        }
        adj.matrix = new_prior_matrix

        ## Perform the likelihood fit if requested
        adj.matrix.fit = lregfit(data,
            adj.matrix,
            adj.matrix.fit,
            regularization,
            command)
    }
    
    
    adj.matrix = list(adj.matrix.pf = adj.matrix.prima.facie,
        adj.matrix.fit = adj.matrix.fit)
    topology = list(adj.matrix = adj.matrix)
    return(topology)
}

# compute the mutual information score for the prim algorithm
compute.mi.score = function ( p_i_j, p_i, p_j ) {
	
	# compute the needed measures
	p_i_not_j = p_i - p_i_j
	p_not_i_j = p_j - p_i_j
	p_not_i_not_j = 1 - p_i - p_j + p_i_j
	
	# compute the 4 terms of mutual information
	m_i_j = p_i_j * log(p_i_j/(p_i*p_j))
	m_i_not_j = p_i_not_j * log(p_i_not_j/(p_i*(1-p_j)))
	m_not_i_j = p_not_i_j * log(p_not_i_j/((1-p_i)*p_j))
	m_not_i_not_j = p_not_i_not_j * log(p_not_i_not_j/((1-p_i)*(1-p_j)))
	
	# NOTE: in our case p_i and p_j are positive numbers in (0,1) with intervals excluded, 
	# hence, the denominator of the scores cannot be 0. 
	# the numerators are numbers in [0,1). So the fraction is a number between [0,Inf). 
	# specifically the log of the fraction is a number in [-Inf,Inf) with Inf excluded and 
	# -Inf when any numerator is 0. 
	# BUT: when numerator is 0, the multiplier of the logarithm is 0 hence givine 0*(-Inf) = NA
	# SO: we will replace any NA with 0 (multiplier wins over logarithm)
	if(is.na(m_i_j)) {
		m_i_j = 0
	}
	if(is.na(m_i_not_j)) {
		m_i_not_j = 0
	}
	if(is.na(m_not_i_j)) {
		m_not_i_j = 0
	}
	if(is.na(m_not_i_not_j)) {
		m_not_i_not_j = 0
	}
	
	# compute the complete mutual information
	mutual_information = m_i_j + m_i_not_j + m_not_i_j + m_not_i_not_j
	
	return(mutual_information)
}


#### end of file -- prim.algorithm.R
