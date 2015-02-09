
library(igraph)
library('Rgraphviz')

t_formula_expansion <- function(input_matrix, atomic_nodes, map) {
  # get node list
  node_list <- colnames(input_matrix)
  # cut input matrix
  min_matrix = input_matrix[-(atomic_nodes+1):-length(node_list), -(atomic_nodes+1):-length(node_list)]
  # create graph from matrix
  min_graph = graph.adjacency(min_matrix)
  # qui ci va il loop!!!
  for (h in ls(map)) {
    min_graph = graph.union(min_graph, map[[h]])
    
    # edge da ricollegare
    h_edge <- input_matrix[h,]
    final_node <- names(h_edge)[which(h_edge==1)]
    
    # nome del nodo
    h_mat <- rowSums(get.adjacency(map[[h]], sparse=FALSE))
    initial_node <- names(h_mat)[which(h_mat==0)]
    
    # rimonta edge mancanti
    for (node in final_node) {
      min_graph <- min_graph + edge(initial_node, node)
    }
    
  }
  min_graph
}

test_t_formula_expansion <- function(){
  h1 <- graph.formula("A"-+"XOR_AB", "B"-+"XOR_AB")
  h2 <- graph.formula("A"-+"AND_AB", "B"-+"AND_AB")
  map <- {}
  map[['h1']] <- h1
  map[['h2']] <- h2
  
  node_list = c("A", "B", "C", "D", "E", "h1", "h2")
  dim_m = length(node_list)
  mat = matrix(0, dim_m, dim_m)
  rownames(mat) = node_list
  colnames(mat) = node_list
  mat["A","C"] = 1
  mat["B", "E"] = 1
  mat["C", "E"] = 1
  mat["D", "E"] = 1
  mat["h1", "C"] = 1
  mat["h2", "D"] = 1
  print(mat)
  print("h1: A XOR B")
  print("h1: A AND B")
  num_atomics = 5
  grafone <- convertGraph(mat, num_atomics, map)
  plot(igraph.to.graphNEL(grafone))
}

