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
  #print(input_matrix)

  # cut input matrix
  margin = length(node_list) - atomic_nodes
  if (length(map) == 0) {
    min_matrix = input_matrix
    return(min_matrix)
  } else {
    min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]
  }
  #print(min_matrix)

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
    #print(hypo)

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
  
  # foreach AND
  min_matrix = get.adjacency(min_graph, sparse = F)
  print(min_matrix)
  and_matrix = NULL
  to_reconnect = list()
  logical_op = list("AND", "OR", "NOT", "XOR")
  for (col in colnames(min_matrix)) {
    # if a node start with AND_.. is already a AND
    prefix = gsub("_.*$", "", col)
    if ( !(prefix %in% logical_op) && sum(min_matrix[,col]) > 1 ) {
      # if colsum > 1 there is a hidden AND and something has to be done..
      print("AND!")
      to_reconnect = append(to_reconnect, col)
      # append a colum to the matrix..
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      # and give her a new name
      colnames(and_matrix)[pos] = paste("AND", col, sep="_")
      
      and_matrix = cbind(and_matrix, matrix(0,nrow = nrow(and_matrix), ncol = 1))
      pos = ncol(and_matrix)
      colnames(and_matrix)[pos] = col
    } else {
      # ..else add a row ad set the correct colname
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      colnames(and_matrix)[pos] = col
    }
  }

  for(row in to_reconnect) {
    and_matrix = rbind(and_matrix, matrix(0, ncol=ncol(and_matrix), nrow = 1))
    pos = nrow(and_matrix)
    rownames(and_matrix)[pos] = paste0("AND_", row)
    and_matrix[paste0("AND_", row),row] = 1
  }
  
  print(and_matrix)
  
  print(colnames(and_matrix))
  print(rownames(and_matrix))
  
  return(and_matrix)
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
