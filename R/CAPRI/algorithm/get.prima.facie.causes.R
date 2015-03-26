#### get.prima.facie.causes.R
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
function(adj.matrix, hypotheses, marginal.probs.distributions, prima.facie.model.distributions, prima.facie.null.distributions, pvalue) {
	
    #structure to save the confidence of the edges
    edge.confidence.matrix <- array(list(), c(2,1));
    edge.confidence.matrix[[1,1]] = array(1, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[2,1]] = array(0, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    
    #verify Suppes' conditions for prima facie causes
    #i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    #verify the temporal priority condition
    temporal.priority = verify.temporal.priority.do.boot(marginal.probs.distributions,pvalue,adj.matrix,edge.confidence.matrix);
    
    #verify the probability raising and background context conditions
    probability.raising = verify.probability.raising.do.boot(prima.facie.model.distributions,prima.facie.null.distributions,pvalue,temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);
    
    #patterns related the the hypotheses
	if(!is.na(hypotheses)[1]) {
		data = list();
		data$hypotheses = hypotheses;
		print(as.patterns(data))
		print(pattern.events(data,as.patterns(data)[1]))
	}
    
    #remove any cycle
    if(length(temporal.priority$not.ordered)>0) {
        weights.matrix = probability.raising$edge.confidence.matrix[[1,1]]+probability.raising$edge.confidence.matrix[[2,1]];
        acyclic.topology = remove.cycles(probability.raising$adj.matrix,weights.matrix,temporal.priority$not.ordered);
        adj.matrix = acyclic.topology$adj.matrix;
    }
    else {
        adj.matrix = probability.raising$adj.matrix;
    }
    
    #save the results and return them
    prima.facie.topology <- list(adj.matrix=adj.matrix,edge.confidence.matrix=probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);

}


#select the best set of prima facie causes per node without bootstrap
#INPUT:
#adj.matrix: adjacency matrix of the initially valid edges
#marginal.probs: marginal probabilities
#prima.facie.model: prima facie model
#prima.facie.null: prima facie null
#RETURN:
#prima.facie.topology: adjacency matrix of the prima facie causes
"get.prima.facie.causes.no.boot" <-
function(adj.matrix, hypotheses, marginal.probs, prima.facie.model, prima.facie.null) {
	
    #verify Suppes' conditions for prima facie causes
    #i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    #verify the temporal priority condition
    temporal.priority = verify.temporal.priority.no.boot(marginal.probs,adj.matrix);
    
    #verify the probability raising and background context conditions
    probability.raising = verify.probability.raising.no.boot(prima.facie.model,prima.facie.null,temporal.priority);
    
    #patterns related the the hypotheses
	if(!is.na(hypotheses[1])) {
		data = list();
		data$hypotheses = hypotheses;
		print(as.patterns(data))
		print(pattern.events(data,as.patterns(data)[1]))
	}
    
    #save the results and return them
    prima.facie.topology = probability.raising;
    return(prima.facie.topology);

}

#### end of file -- get.prima.facie.causes.R
