#### verify.temporal.priority.do.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#verify the temporal priority condition
#INPUT:
#marginal.probs.distributions: distributions of the bootstrapped marginal probabilities
#pvalue: minimum pvalue for the Mann-Whitney U tests to be significant
#adj.matrix: adjacency matrix of the topology
#edge.confidence.matrix: matrix of the confidence of each edge
#RETURN:
#temporal.priority: list describing the causes where temporal priority is verified
"verify.temporal.priority.do.boot" <-
function(marginal.probs.distributions, pvalue, adj.matrix, edge.confidence.matrix) {
    #evalutate the temporal priority condition for each par of edges
    not.ordered = list();
    counter = 0;
    for(i in 1:nrow(adj.matrix)) {
        for(j in i:ncol(adj.matrix)) {
            #the diagonal (self cause) and the other invalid edges have not to be considered
            if(adj.matrix[i,j]!=0 || adj.matrix[j,i]!=0) {
                #[i,j] refers to causation i --> j
                #temporal priority condition: if P(i)>P(j) the edge i --> j is valid for temporal priority
                #test i --> j
                first.pvalue.i.j = wilcox.test(unlist(marginal.probs.distributions[i,1]),unlist(marginal.probs.distributions[j,1]),alternative="greater",mu=0)$p.value;
                if(is.na(first.pvalue.i.j) || is.nan(first.pvalue.i.j)) {
                    #in this case the two distributions are exactly identical
                    first.pvalue.i.j = 1;
                }
                #test j --> i
                first.pvalue.j.i = wilcox.test(unlist(marginal.probs.distributions[j,1]),unlist(marginal.probs.distributions[i,1]),alternative="greater",mu=0)$p.value;
                if(is.na(first.pvalue.j.i) || is.nan(first.pvalue.j.i)) {
                    #in this case the two distributions are exactly identical
                    first.pvalue.j.i = 1;
                }
                #in this case i is before j and j --> i is not valid
                if(first.pvalue.j.i>=pvalue && first.pvalue.i.j<pvalue) {
                    #[j,i] = 0 means j is after i, i.e. it can not be causing i
                    adj.matrix[j,i] = 0;
                }
                #in this case j is before i and i --> j is not valid
                else if(first.pvalue.j.i<pvalue && first.pvalue.i.j>=pvalue) {
                    #[i,j] = 0 means i is after j, i.e. it can not be causing j
                    adj.matrix[i,j] = 0;
                }
                #in this case, a total time order between i and j can not be defined
                else {
                    #no temporal priority induced by the topology can be inferred
                    counter = counter + 1;
                    curr.not.ordered = array(-1, c(2,1));
                    curr.not.ordered[1,1] = i;
                    curr.not.ordered[2,1] = j;
                    not.ordered[counter] = list(curr.not.ordered);
                }
                #save the confidence for i --> j and j --> i
                tmp = edge.confidence.matrix[[1,1]];
                tmp[i,j] = first.pvalue.i.j;
                tmp[j,i] = first.pvalue.j.i;
                edge.confidence.matrix[1,1] = list(tmp);
            }
            else {
                tmp = edge.confidence.matrix[[1,1]];
                tmp[i,j] = 1;
                tmp[j,i] = 1;
                edge.confidence.matrix[1,1] = list(tmp);
            }
        }
    }
    temporal.priority <- list(adj.matrix=adj.matrix,edge.confidence.matrix=edge.confidence.matrix,not.ordered=not.ordered);
    return(temporal.priority);
}

#### end of file -- verify.temporal.priority.do.boot.R
