#### verify.probability.raising.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# verify the probability raising condition
# INPUT:
# prima.facie.model.distributions: distributions of the prima facie model
# prima.facie.null.distributions: distributions of the prima facie null
# pvalue: minimum pvalue for the Mann-Whitney U tests to be significant
# adj.matrix: adjacency matrix of the topology
# edge.confidence.matrix: matrix of the confidence of each edge
# RETURN:
# probability.raising: list describing the causes where probability raising is verified
"verify.probability.raising.do.boot" <-
function( prima.facie.model.distributions, prima.facie.null.distributions, pvalue, adj.matrix, edge.confidence.matrix ) {
	
    # evaluate the probability raising condition
    for(i in 1:nrow(adj.matrix)) {
        for(j in i:ncol(adj.matrix)) {
        	
            # the diagonal (self cause) and the other invalid edges have not to be considered
            if(adj.matrix[i,j]!=0 || adj.matrix[j,i]!=0) {
            	
                # pvalue for the probability raising condition for i --> j
                second.pvalue.i.j = wilcox.test(unlist(prima.facie.model.distributions[i,j]),unlist(prima.facie.null.distributions[i,j]),alternative="greater",mu=0)$p.value;
                if(is.na(second.pvalue.i.j) || is.nan(second.pvalue.i.j)) {
                    # in this case the two distributions are exactly identical
                    second.pvalue.i.j = 1;
                }
                
                # in this case i --> j is not valid
                if(second.pvalue.i.j>=pvalue) {
                    adj.matrix[i,j] = 0;
                }
                
                # pvalue for the probability raising condition for j --> i
                second.pvalue.j.i = wilcox.test(unlist(prima.facie.model.distributions[j,i]),unlist(prima.facie.null.distributions[j,i]),alternative="greater",mu=0)$p.value;
                if(is.na(second.pvalue.j.i) || is.nan(second.pvalue.j.i)) {
                    # in this case the two distributions are exactly identical
                    second.pvalue.j.i = 1;
                }
                
                # in this case j --> i is not valid
                if(second.pvalue.j.i>=pvalue) {
                    adj.matrix[j,i] = 0;
                }
                
                # save the confidence for i-->j and j --> i
                tmp = edge.confidence.matrix[[2,1]];
                tmp[i,j] = second.pvalue.i.j;
                tmp[j,i] = second.pvalue.j.i;
                edge.confidence.matrix[2,1] = list(tmp);
                
            }
            else {
                tmp = edge.confidence.matrix[[2,1]];
                tmp[i,j] = 1;
                tmp[j,i] = 1;
                edge.confidence.matrix[2,1] = list(tmp);
            }
            
        }
        
    }
    
    # save the results and return them
    probability.raising <- list(adj.matrix=adj.matrix,edge.confidence.matrix=edge.confidence.matrix);
    return(probability.raising);

}


# verify the probability raising condition without bootstrap
# INPUT:
# prima.facie.model: prima facie model
# prima.facie.null: prima facie null
# adj.matrix: adjacency matrix of the topology
# edge.confidence.matrix: matrix of the confidence of each edge
# RETURN:
# probability.raising: adjacency matrix where temporal priority is verified
"verify.probability.raising.no.boot" <-
function( prima.facie.model, prima.facie.null, adj.matrix, edge.confidence.matrix ) {
	
    # evaluate the probability raising condition
    for(i in 1:nrow(adj.matrix)) {
        for(j in i:ncol(adj.matrix)) {
        	
            # the diagonal (self cause) and the other invalid edges have not to be considered
            # probability raising condition: if P(j|i)>P(j|not i) the edge i --> j is valid for probability raising
            if(adj.matrix[i,j]!=0 || adj.matrix[j,i]!=0) {
            	
                # in this case i --> j is not valid
                if(prima.facie.model[i,j]<=prima.facie.null[i,j]) {
                    adj.matrix[i,j] = 0;
                }
                
                # in this case j --> i is not valid
                if(prima.facie.model[j,i]<=prima.facie.null[j,i]) {
                    adj.matrix[j,i] = 0;
                }
                
                # save the confidence for i-->j and j --> i
                tmp = edge.confidence.matrix[[2,1]];
                tmp[i,j] = min(prima.facie.null[i,j]/prima.facie.model[i,j],1);
                tmp[j,i] = min(prima.facie.null[j,i]/prima.facie.model[j,i],1);
                edge.confidence.matrix[2,1] = list(tmp);
                
            }
            else {
                tmp = edge.confidence.matrix[[2,1]];
                tmp[i,j] = 1;
                tmp[j,i] = 1;
                edge.confidence.matrix[2,1] = list(tmp);
            }
            
        }
    }
    
    # save the results and return them
    probability.raising <- list(adj.matrix=adj.matrix,edge.confidence.matrix=edge.confidence.matrix);
    return(probability.raising);
}

#### end of file -- verify.probability.raising.R
