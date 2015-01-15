#### verify.temporal.priority.no.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#verify the temporal priority condition without bootstrap
#INPUT:
#marginal.probs: marginal probabilities
#adj.matrix: adjacency matrix of the topology
#RETURN:
#temporal.priority: adjacency matrix where temporal priority is verified
"verify.temporal.priority.no.boot" <-
function(marginal.probs, adj.matrix) {
    for(i in 1:nrow(adj.matrix)) {
        for(j in i:ncol(adj.matrix)) {
            #the diagonal (self cause) and the other invalid edges have not to be considered
            if(adj.matrix[i,j]!=0 || adj.matrix[j,i]!=0) {
                #[i,j] refers to causation i --> j
                #temporal priority condition: if P(i)>P(j) the edge i --> j is valid for temporal priority
                if(marginal.probs[i,1]>marginal.probs[j,1]) {
                    #in this case j --> i is not valid
                    adj.matrix[j,i] = 0;
                }
                else {
                    #in this case i --> j is not valid
                    adj.matrix[i,j] = 0;
                }
            }
        }
    }
    temporal.priority = adj.matrix;
    return(temporal.priority);
}

#### end of file -- verify.temporal.priority.no.boot.R
