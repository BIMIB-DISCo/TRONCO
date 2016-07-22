#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

# reconstruct the best topology based on probabilistic causation and maximum likelihood estimation
# @title gabow.fit
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param score the score to be used, could be either pointwise mutual information (pmi) or conditional entropy (entropy)
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent should I be verbose?
# @param epos error rate of false positive errors
# @param eneg error rate of false negative errors
# @return topology: the reconstructed tree topology
#
gabow.fit <- function(dataset,
                      regularization = "no_reg",
                      score = "pmi",
                      do.boot = TRUE,
                      nboot = 100,
                      pvalue = 0.05,
                      min.boot = 3,
                      min.stat = TRUE,
                      boot.seed = NULL,
                      silent = FALSE,
                      epos = 0.0,
                      eneg = 0.0,
                      do.raising = TRUE ) {

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
    
    valid.dataset = check.dataset(dataset, adj.matrix, FALSE, epos, eneg)
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
                                            silent,
                                            epos,
                                            eneg);
    } else {
        if (!silent)
            cat('*** Computing selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.no.boot(dataset,
                                            NA,
                                            adj.matrix,
                                            silent,
                                            epos,
                                            eneg);
    }

    ## Add back in any connection invalid for the probability raising
    ## theory.
    
    if (length(invalid.events) > 0) {
        # save the correct acyclic matrix
        adj.matrix.cyclic.tp.valid = prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp
        adj.matrix.cyclic.valid = prima.facie.parents$adj.matrix$adj.matrix.cyclic
        adj.matrix.acyclic.valid = prima.facie.parents$adj.matrix$adj.matrix.acyclic
        for (i in 1:nrow(invalid.events)) {
            prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1
            prima.facie.parents$adj.matrix$adj.matrix.cyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1
            prima.facie.parents$adj.matrix$adj.matrix.acyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1
        }
        # if the new cyclic.tp contains cycles use the previously computed matrix
        if (!is.dag(graph.adjacency(prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp))) {
            prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp = adj.matrix.cyclic.tp.valid
        }
        # if the new cyclic contains cycles use the previously computed matrix
        if (!is.dag(graph.adjacency(prima.facie.parents$adj.matrix$adj.matrix.cyclic))) {
            prima.facie.parents$adj.matrix$adj.matrix.cyclic = adj.matrix.cyclic.valid
        }
        # if the new acyclic contains cycles use the previously computed matrix
        if (!is.dag(graph.adjacency(prima.facie.parents$adj.matrix$adj.matrix.acyclic))) {
            prima.facie.parents$adj.matrix$adj.matrix.acyclic = adj.matrix.acyclic.valid
        }
    }
    
    adj.matrix.prima.facie =
        prima.facie.parents$adj.matrix$adj.matrix.acyclic
        
    # set the input prior matrix for Gabow
    if(do.raising==TRUE) {
        adj.matrix.gabow = prima.facie.parents$adj.matrix$adj.matrix.cyclic
    }
    else {
        adj.matrix.gabow = prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp
    }

    ## Perform the likelihood fit with the required strategy.
    
    model = list();
    for(my_score in score) {
        
        # as Gabow algorithm is computationally expensive,
        # I first reconstruct with "no_reg" for each score
        best.parents =
                perform.likelihood.fit.gabow(dataset,
                                       adj.matrix.gabow,
                                       regularization = "no_reg",
                                       score = my_score,
                                       marginal.probs = prima.facie.parents$marginal.probs,
                                       joint.probs = prima.facie.parents$joint.probs)
        
        # if Gabow was not performed due to memory issues,
        # use Edmonds instead
        if(is.null(best.parents)) {
            best.parents =
                perform.likelihood.fit.edmonds(dataset,
                                       adj.matrix.prima.facie,
                                       regularization = reg,
                                       score = my_score,
                                       marginal.probs = prima.facie.parents$marginal.probs,
                                       joint.probs = prima.facie.parents$joint.probs)
        }
                                       
        # now I perform the likelihood fit with each regularizator
        for (reg in regularization) {

            ## Perform the likelihood fit with the chosen regularization
            ## score on the prima facie topology.
            
            if (!silent)
                cat('*** Performing likelihood-fit with regularization:', reg, '.\n')
                
            
            # adjacency matrix of the topology reconstructed by likelihood fit
            curr_adj.matrix.fit = array(0,c(nrow(best.parents$adj.matrix$adj.matrix.fit),ncol(best.parents$adj.matrix$adj.matrix.fit)))
            rownames(curr_adj.matrix.fit) = colnames(dataset)
            colnames(curr_adj.matrix.fit) = colnames(dataset)
            
            curr_adj.matrix.fit = lregfit(as.categorical.dataset(dataset),
                                    best.parents$adj.matrix$adj.matrix.fit,
                                    curr_adj.matrix.fit,
                                    reg,
                                    command = "hc")
                                    
                                    
            ## Save the results and return them.
            adj.matrix_reg = list(adj.matrix.pf = adj.matrix.gabow,
                                  adj.matrix.fit = curr_adj.matrix.fit)
            curr_best.parents = list(adj.matrix = adj.matrix_reg)
            
            ## Set the structure to save the conditional probabilities of
            ## the reconstructed topology.
    
            reconstructed.model = create.model(dataset,
                curr_best.parents,
                prima.facie.parents)
    
            model.name = paste('gabow', reg, my_score, sep='_')
            model[[model.name]] = reconstructed.model
            
        }
        
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "GABOW",
             regularization = regularization,
             score = score,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent,
             error.rates = list(epos=epos,eneg=eneg),
             do.raising = do.raising)

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


# reconstruct the best causal topology by Gabow algorithm combined with probabilistic causation
# @title perform.likelihood.fit.gabow
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param regularization regularization term to be used in the likelihood fit
# @param command type of search, either hill climbing (hc) or tabu (tabu)
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.gabow = function(dataset,
                                        adj.matrix,
                                        regularization,
                                        score,
                                        command = "hc",
                                        marginal.probs,
                                        joint.probs) {

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
    
    # detect the acyclic parts of the graph,
    # being the strongly connected components of degree 1
    my_graph = graph_from_adjacency_matrix(adj.matrix)
    strongly_connected = clusters(my_graph,mode="strong")

    # estimate the size required to save the permutations for the 
    # bigger strongly connected component
    max_component = max(strongly_connected$csize)
    max_rows = factorial(max_component)
    max_columns = max_component
    max_num_entries = max_rows * max_columns
    max_memory_size = (max_num_entries)*4/1024^3
    
    # perform Gabow only if the maximum strongly connected component 
    # requires allocation of at most 10 GBs of memory
    if(max_memory_size<=5) {
        valid_nodes = c()
        for(i in 1:strongly_connected$no) {
            if(length(which(strongly_connected$membership==i))==1) {
                valid_nodes = c(valid_nodes,which(strongly_connected$membership==i))
            }
        }
           
        # set at most one parent per node based on mutual information
        # only for the acyclic parts of the graph
        if(length(valid_nodes)>0) {
            for (i in valid_nodes) {
                    
                # consider the parents of i
                curr_parents = which(adj.matrix[,i] == 1)
                
                # if I have more then one valid parent
                if (length(curr_parents) > 1) {
                    
                    # find the best parent
                    curr_best_parent = -1
                    curr_best_score = NA
                    for (j in curr_parents) {
                        
                        # if the event is valid
                        if(joint.probs[i,j]>=0) {
                            # compute the chosen score
                            new_score = compute.edmonds.score(joint.probs[i,j],marginal.probs[i],marginal.probs[j],score)
                        }
                        # else, if the two events are indistinguishable
                        else if(joint.probs[i,j]<0) {
                            new_score = Inf
                        }
                        
                        if (is.na(curr_best_score) || new_score > curr_best_score) {
                            curr_best_parent = j
                            curr_best_score = new_score
                        }
                    }
                    
                    # set the best parent
                    for (j in curr_parents) {
                        if (j != curr_best_parent) {
                            adj.matrix[j,i] = 0
                        }
                    }
                }
            }
        }
           
        # now consider the acyclic parts of the graph
        if(length(valid_nodes)<nrow(adj.matrix)) {
            # now I consider each strongly connected componet with size greater then 1
            for(i in 1:strongly_connected$no) {
                if(length(which(strongly_connected$membership==i))>1) {
                    adj.matrix = find.best.subtree(adj.matrix,which(strongly_connected$membership==i),marginal.probs,joint.probs,score)
                }
            }
        }
        
        # perform the likelihood fit if requested
        adj.matrix.fit = lregfit(data,
            adj.matrix,
            adj.matrix.fit,
            regularization,
            command)
        
        ## Save the results and return them.
        
        adj.matrix =
            list(adj.matrix.pf = adj.matrix.prima.facie,
                 adj.matrix.fit = adj.matrix.fit)
        topology = list(adj.matrix = adj.matrix)
    }
    else {
        warning("Too big strongly connected components: Gabow can not be performed, using Edmonds instead.")
        topology = NULL
    }
        
    return(topology)

}

find.best.subtree = function( adj.matrix, subtree.nodes, marginal.probs, joint.probs, score ) {
    
    # enumerate the possible orderings within subtree.nodes
    possible.orderings = permutations(length(subtree.nodes),length(subtree.nodes))
    
    # consider all the possible trees within subtree.nodes
    all_tree_scores = sapply(1:nrow(possible.orderings),FUN=function(x) {
        curr.ordering = subtree.nodes[possible.orderings[x,]]
        curr.adj.matrix = adj.matrix
        # remove any loop accordingly to this ordering
        for(j in 2:length(curr.ordering)) {
            curr_invalid = curr.ordering[1:(j-1)]
            curr.adj.matrix[curr.ordering[j],curr_invalid] = 0
        }
        # compute the score of curr.adj.matrix
        curr_score = get.best.scored.tree(curr.adj.matrix,subtree.nodes,marginal.probs,joint.probs,score)$score
        return(curr_score)
    })
    
    # pick the best scored tree
    best.ordering = subtree.nodes[possible.orderings[which(all_tree_scores==max(all_tree_scores))[1],]]
    
    # remove any loop accordingly to this ordering
    best_adj.matrix = adj.matrix
    for(j in 2:length(best.ordering)) {
        curr_invalid = best.ordering[1:(j-1)]
        best_adj.matrix[best.ordering[j],curr_invalid] = 0
    }
    
    # compute the score of curr.adj.matrix
    best_adj.matrix = get.best.scored.tree(best_adj.matrix,subtree.nodes,marginal.probs,joint.probs,score)$adj.matrix
    
    return(best_adj.matrix)
}

get.best.scored.tree = function( adj.matrix, subtree.nodes, marginal.probs, joint.probs, score_type ) {
    
    score = 0
    total.arcs = 0
    attached.arcs = 0
    
    # set at most one parent per node based on mutual information
    for (i in subtree.nodes) {
        
        # consider the parents of i
        curr_parents = which(adj.matrix[,i] == 1)
        
        # if I have more then one valid parent
        if (length(curr_parents) > 1) {
            
            # find the best parent
            curr_best_parent = -1
            curr_best_score = NA
            for (j in curr_parents) {
                
                # if the event is valid
                if(joint.probs[i,j]>=0) {
                    # compute the chosen score
                    new_score = compute.edmonds.score(joint.probs[i,j],marginal.probs[i],marginal.probs[j],score_type)
                }
                # else, if the two events are indistinguishable
                else if(joint.probs[i,j]<0) {
                    new_score = Inf
                }
                
                if (is.na(curr_best_score) || new_score > curr_best_score) {
                    curr_best_parent = j
                    curr_best_score = new_score
                }
            }
            
            if(curr_best_score!=Inf) {
                score = score + curr_best_score
                total.arcs = total.arcs + 1
            }
            
            # set the best parent
            invalid.curr.parents = which(curr_parents!=curr_best_parent)
            if(length(invalid.curr.parents)>0) {
                adj.matrix[curr_parents[invalid.curr.parents],i] = 0
            }
            
            attached.arcs = attached.arcs + 1
            
        }
    }
    
    # if I have an empty topology, it has the worst score
    if(attached.arcs==0) {
        score = -Inf
    }
    # if I have only arcs between indistinguishable events
    else if(attached.arcs>0 && total.arcs==0) {
        score = Inf
    }
    # if I have a valid topology which I want to weight
    # for the number of arcs
    else {
        mean_score = score / total.arcs
        score = score + (mean_score * (attached.arcs-total.arcs))
    }
    
    res = list(adj.matrix=adj.matrix,score=score)
    
    return(res)
    
}


#### end of file -- gabow.algorithm.R
