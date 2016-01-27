#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Convert a TRONCO object in a Bnlean network
#' @title as.bnlearn.network
#'
#' @examples
#' data(test_model)
#' as.bnlearn.network(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#'
as.bnlearn.network <- function(data, regularization = "bic") {

	## Check if there is a reconstructed model.

    if(!has.model(data)) {
        stop('This dataset doesn\'t have.')
    }

    ## Check if the selected regularization is used in the model.

    if (!regularization %in% names(data$model)) {
        stop(paste(regularization, "not in model"))
    }

    ## Get genotypes and data.

    genotypes = as.genotypes(data)
    genotypes = keysToNames(data, genotypes)
	names(colnames(genotypes)) = NULL

    adj.matrix = get(regularization, as.adj.matrix(data))
	adj.matrix = keysToNames(data, adj.matrix)
	names(colnames(adj.matrix)) = NULL
	names(rownames(adj.matrix)) = NULL
	
	bayes.net = NULL
	    	
    ## Create a categorical data frame from the dataset.
    df = array("missing",c(nrow(genotypes),ncol(genotypes)))
    for (i in 1:nrow(genotypes)) {
        for (j in 1:ncol(genotypes)) {
            if(genotypes[i,j]==1) {
                df[i,j] = "observed"
            }
        }
    }
    df = as.data.frame(df)
    my.names = names(df)
    
    for (i in 1:length(my.names)) {
        my.names[i] = toString(i)
    }
    colnames(df) = colnames(genotypes)
    bayes.net$data = df
        
    ## Create the Bayesian Network of the fitted model.
    bayes.net$net = empty.graph(colnames(genotypes))
    for (i in 1:nrow(adj.matrix)) {
    	for(j in 1:ncol(adj.matrix)) {
    		if(adj.matrix[i,j]==1) {
    			bayes.net$net = set.arc(
    				bayes.net$net, 
    				from=colnames(genotypes)[i], 
    				to=colnames(genotypes)[j])
    		}
    	}
    }
   	return(bayes.net) 
}
