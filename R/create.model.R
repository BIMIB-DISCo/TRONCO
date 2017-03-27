#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2017, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# create a reconstructed model
# param datataset genotypes matrix
# param best.parents result of perform.likelihood.fit.<alg>
# param prima.facie.parents result of get.prima.facie.parents <boot|no boot>
create.model <- function(dataset,
                         best.parents,
                         prima.facie.parents){

    ## Set the structure to save the conditional probabilities of
    ## the reconstructed topology.
    
    parents.pos.fit = array(list(),c(ncol(dataset),1));
    conditional.probs.fit = array(list(),c(ncol(dataset),1));

    ## Compute the conditional probabilities.
    
    for (i in 1:ncol(dataset)) {
        for (j in 1:ncol(dataset)) {
            if (i!=j && best.parents$adj.matrix$adj.matrix.fit[i, j] == 1) {
                parents.pos.fit[j,1] =
                    list(c(unlist(parents.pos.fit[j, 1]), i))
                conditional.probs.fit[j,1] =
                    list(c(unlist(conditional.probs.fit[j, 1]),
                           prima.facie.parents$joint.probs[i, j] /
                               prima.facie.parents$marginal.probs[i]))
            }
        }
    }
    parents.pos.fit[unlist(lapply(parents.pos.fit, is.null))] = list(-1)
    conditional.probs.fit[unlist(lapply(conditional.probs.fit, is.null))] = list(1)

    ## Perform the estimation of the probabilities if requested.
    
    estimated.error.rates.fit =
        list(error.fp = NA,
             error.fn = NA)
    estimated.probabilities.fit =
        list(marginal.probs = NA,
             joint.probs = NA,
             conditional.probs = NA)

    ## Set results for the current regolarizator
    probabilities.observed =
        list(marginal.probs = prima.facie.parents$marginal.probs,
             joint.probs = prima.facie.parents$joint.probs,
             conditional.probs = conditional.probs.fit)
    probabilities.fit =
        list(estimated.marginal.probs = estimated.probabilities.fit$marginal.probs,
             estimated.joint.probs = estimated.probabilities.fit$joint.probs,
             estimated.conditional.probs = estimated.probabilities.fit$conditional.probs)
    probabilities =
        list(probabilities.observed = probabilities.observed,
             probabilities.fit = probabilities.fit)
    parents.pos = parents.pos.fit
    error.rates = estimated.error.rates.fit

    ## Save the results for the model.

    result =
        list(probabilities = probabilities,
             parents.pos = parents.pos,
             error.rates = error.rates,
             adj.matrix = best.parents$adj.matrix)

    return(result)
}



#### end of file -- create.model.R
