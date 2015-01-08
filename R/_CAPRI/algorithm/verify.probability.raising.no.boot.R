#### verify.probability.raising.no.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#verify the probability raising condition without bootstrap
#INPUT:
#prima.facie.model: prima facie model
#prima.facie.null: prima facie null
#adj.matrix: adjacency matrix of the topology
#RETURN:
#probability.raising: adjacency matrix where temporal priority is verified
"verify.probability.raising.no.boot" <-
function(prima.facie.model, prima.facie.null, adj.matrix) {
    for(i in 1:nrow(adj.matrix)) {
        for(j in i:ncol(adj.matrix)) {
            #the diagonal (self cause) and the other invalid edges have not to be considered
            #temporal priority condition: if P(j|i)>P(j|not i) the edge i --> j is valid for temporal priority
            if(adj.matrix[i,j]!=0) {
                #verify i --> j
                if(prima.facie.model[i,j]<=prima.facie.null[i,j]) {
                    adj.matrix[i,j] = 0;
                }
            }
            if(adj.matrix[j,i]!=0) {
                #verify j --> i
                if(prima.facie.model[j,i]<=prima.facie.null[j,i]) {
                    adj.matrix[j,i] = 0;
                }
            }
        }
    }
    probability.raising = adj.matrix;
    return(probability.raising);
}

#### end of file -- verify.probability.raising.no.boot.R
