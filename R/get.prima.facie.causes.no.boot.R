#### get.prima.facie.causes.no.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select the best set of prima facie causes per node without bootstrap
#INPUT:
#marginal.probs: marginal probabilities
#prima.facie.model: prima facie model
#prima.facie.null: prima facie null
#RETURN:
#prima.facie.topology: adjacency matrix of the prima facie causes
"get.prima.facie.causes.no.boot" <-
function(marginal.probs, prima.facie.model, prima.facie.null) {
    #structure to save the adjacency matrix of the prima facie causes
    #the matrix is complete: all the edges are plausible causations except for the diagonal ones (self cause)
    adj.matrix <- array(1, c(ncol(prima.facie.model),ncol(prima.facie.model)));
    adj.matrix[row(adj.matrix)==col(adj.matrix)] <- 0;
    #verify Suppes' conditions for prima facie causes
    #i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    #verify the temporal priority condition
    temporal.priority = verify.temporal.priority.no.boot(marginal.probs,adj.matrix);
    #verify the probability raising and background context conditions
    probability.raising = verify.probability.raising.no.boot(prima.facie.model,prima.facie.null,temporal.priority);
    prima.facie.topology = probability.raising;
    return(prima.facie.topology);
}

#### end of file -- get.prima.facie.causes.no.boot.R
