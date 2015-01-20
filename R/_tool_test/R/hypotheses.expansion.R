#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

hypotheses.expansion <- function(input_matrix, 
                                atomic_nodes = NULL, 
                                map = list()) {
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  
  # get node list
  node_list <- colnames(input_matrix)
  print(input_matrix)

  # cut input matrix
  margin = length(node_list) - atomic_nodes
  if (length(map) == 0) {
    min_matrix = input_matrix
    return(min_matrix)
  } else {
    min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]
  }
  print(min_matrix)

  # create graph from matrix
  min_graph = graph.adjacency(min_matrix)
  

  # foreach hypothesis
  for (h in ls(map)) {
    if(length(which(input_matrix[h,] == 1)) == 0) {
      break
    }
    
    # eros! please give me the transposed matrix
    hypo = map[[h]]
    
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

hypo.plot = function(capri, data, hypotheses = NULL, font=14) {
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    install.packages('Rgraphviz', dependencies = TRUE)
    library(Rgraphviz)
  }
  
  c_matrix = capri$adj.matrix$adj.matrix.bic
  
  colnames(c_matrix) = colnames(capri$dataset);
  rownames(c_matrix) = colnames(capri$dataset);
  
  if (is.null(hypotheses)) {
    hstruct = NULL
    num_h = NULL
  } else {
    hstruct = hypotheses$hypotheses$hstructure
    num_h = length(hstruct)
  }

  hypo_mat = hypotheses.expansion(c_matrix, num_h, hstruct)
  hypo_graph = graph.adjacency(hypo_mat)
  v_names = gsub("_.*$", "", V(hypo_graph)$name)
  new_name = list()
  for(v in v_names) {
    if(v %in% rownames(data$annotations)) {
      n = data$annotations[v,"event"]
      t = data$annotations[v,"type"]
      new_name = append(new_name, paste(n, t))
    } else {
      new_name = append(new_name, v)
    }
  }
  V(hypo_graph)$label = new_name
  graph <- igraph.to.graphNEL(hypo_graph)
  z = V(hypo_graph)$label
  names(z) = nodes(graph)
  nAttrs = list()
  nAttrs$label = z

  # nAttrs$fontsize = rep('8', length(nAttrs$label))
  # names(nAttrs$fontsize) = z

  attrs <- list(node = list(fixedsize = FALSE, fontsize=font)) 
  # attrs$node$fontsize=8 
   
  # print(nAttrs)
   
  plot(graph, nodeAttrs=nAttrs, attrs=attrs)
}
