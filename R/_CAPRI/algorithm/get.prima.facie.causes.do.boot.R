#### get.prima.facie.causes.do.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select the best set of prima facie causes per node
#INPUT:
#adj.matrix: adjacency matrix of the initially valid edges
#marginal.probs.distributions: distributions of the bootstrapped marginal probabilities
#prima.facie.model.distributions: distributions of the prima facie model
#prima.facie.null.distributions: distributions of the prima facie null
#pvalue: minimum pvalue for the Mann-Whitney U tests to be significant
#RETURN:
#prima.facie.topology: list describing the topology of the prima facie causes
"get.prima.facie.causes.do.boot" <-
function(adj.matrix, marginal.probs.distributions, prima.facie.model.distributions, prima.facie.null.distributions, pvalue) {
    #structure to save the confidence of the edges
    edge.confidence.matrix <- array(list(), c(2,1));
    edge.confidence.matrix[[1,1]] = array(1, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[2,1]] = array(0, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    #verify Suppes' conditions for prima facie causes
    #i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    #verify the temporal priority condition
    temporal.priority = verify.temporal.priority(marginal.probs.distributions,pvalue,adj.matrix,edge.confidence.matrix);
    #verify the probability raising and background context conditions
    probability.raising = verify.probability.raising(prima.facie.model.distributions,prima.facie.null.distributions,pvalue,temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);
    #remove any cycle
    if(length(temporal.priority$not.ordered)>0) {
        weights.matrix = probability.raising$edge.confidence.matrix[[1,1]]+probability.raising$edge.confidence.matrix[[2,1]];
        acyclic.topology = remove.cycles(probability.raising$adj.matrix,weights.matrix,temporal.priority$not.ordered);
        adj.matrix = acyclic.topology$adj.matrix;
    }
    else {
        adj.matrix = probability.raising$adj.matrix;
    }
    prima.facie.topology <- list(adj.matrix=adj.matrix,edge.confidence.matrix=probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);
}

#### end of file -- get.prima.facie.causes.do.boot.R
 