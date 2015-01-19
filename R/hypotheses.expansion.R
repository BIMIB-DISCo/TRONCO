#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

library('Rgraphviz')

hypotheses.expansion <- function(input_matrix, 
                                atomic_nodes, 
                                map) {
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  
  # get node list
  node_list <- colnames(input_matrix)
  print(input_matrix)

  # cut input matrix
  margin = length(node_list) - atomic_nodes
  min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]
  print(min_matrix)

  # create graph from matrix
  min_graph = graph.adjacency(min_matrix)
  

  # foreach hypothesis
  for (h in ls(map)) {
    if(length(which(input_matrix[h,] == 1)) == 0) {
      break
    }
    
    # eros! please give me the transposed matrix
    hypo = t(map[[h]])
    
    # create graph from hypo
    hypo_graph = graph.adjacency(hypo)
    print(hypo)

    # add this graph to main graph
    min_graph = graph.union(min_graph, hypo_graph)
    
    # edge to reconstruct
    h_edge <- input_matrix[h,]
    final_node <- names(h_edge)[which(h_edge==1)]
    
    # name of this node
    h_mat <- rowSums(get.adjacency(hypo_graph, sparse=FALSE))
    initial_node <- names(h_mat)[which(h_mat==0)]
    
    # recreate lost edge
    for (node in final_node) {
      min_graph <- min_graph + edge(initial_node, node)
    }
    
  }
  return(get.adjacency(min_graph, sparse = F))
}
