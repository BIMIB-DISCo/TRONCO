#### capri.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

hypotheses.expansion <- function(input_matrix, 
                                 map = list()) {
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  num_hypos = length(map)
  
  # get node list
  node_list <- colnames(input_matrix)
  #print(input_matrix)
  
  # cut input matrix
  margin = length(node_list) - num_hypos
  
  # check if there are hypotheses
  if (num_hypos == 0) {
    # if no hypos do nothings..
    min_matrix = input_matrix
  } else {
    # ..else expand them
    min_matrix = input_matrix[-(margin+1):-length(node_list), -(margin+1):-length(node_list)]
    
    # create graph from matrix
    min_graph = graph.adjacency(min_matrix)
    
    
    # foreach hypothesis
    # print(ls(map))
    # for (h in ls(map)) {
    # print(h)
    # }
    
    
    for (h in ls(map)) {
      
      # print(input_matrix[h,])
      
      # print(any(input_matrix[h,] == 1))
      
      
      
      if(length(which(input_matrix[h,] == 1)) == 0) {
        next
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
    min_matrix = get.adjacency(min_graph, sparse = F)
  }
  
  # now expand the hidden AND
  # print(min_matrix)
  
  and_matrix = NULL
  to_reconnect = list()
  logical_op = list("AND", "OR", "NOT", "XOR")
  
  # foreach AND column
  for (col in colnames(min_matrix)) {
    prefix = gsub("_.*$", "", col)
    if ( !(prefix %in% logical_op) && sum(min_matrix[,col]) > 1 ) {
      # not logical operator and colsum > 1 there is a hidden AND and something has to be done..
      
      # remember to reconnect the fake and to this node
      to_reconnect = append(to_reconnect, col)
      # append a colum from the old matrix..
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      # and give her a new name based on the old one
      colnames(and_matrix)[pos] = paste("AND", col, sep="_")
      
      # append a 0 columl to the matrix..
      and_matrix = cbind(and_matrix, matrix(0,nrow = nrow(and_matrix), ncol = 1))
      pos = ncol(and_matrix)
      # and give her the old name
      colnames(and_matrix)[pos] = col
    } else {
      # ..else add the row taken from the old matrix and set the correct colname
      and_matrix = cbind(and_matrix, min_matrix[,col])
      pos = ncol(and_matrix)
      colnames(and_matrix)[pos] = col
    }
  }
  
  # now reconnect AND node to his gene (AND_Gene7 -> Gene7)
  for(row in to_reconnect) {
    and_matrix = rbind(and_matrix, matrix(0, ncol=ncol(and_matrix), nrow = 1))
    pos = nrow(and_matrix)
    rownames(and_matrix)[pos] = paste0("AND_", row)
    and_matrix[paste0("AND_", row),row] = 1
  }
  
  # sort col and row (igraph wants the same order)
  and_matrix = and_matrix[,order(colnames(and_matrix))]
  and_matrix = and_matrix[order(rownames(and_matrix)),]
  
  # print(and_matrix)  
  return(and_matrix)
}


hypo.plot = function(data, font=14, pf = FALSE, disconnected=FALSE, name=deparse(substitute(capri))) {
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    install.packages('Rgraphviz', dependencies = TRUE)
    library(Rgraphviz)
  }
  
  c_matrix = data$adj.matrix$adj.matrix.bic
  
  
  
  if(pf) c_matrix = data$adj.matrix$adj.matrix.prima.facie
  
  # the TRONCO object
  data = data$data
  
  # hypotheses
  hypotheses = data$data$hypotheses
  
  if (is.null(hypotheses)) {
    hstruct = NULL
    num_h = NULL
  } else {
    hstruct = hypotheses$hstructure
    num_h = length(hstruct)
  }
  
  hypo_mat = hypotheses.expansion(c_matrix, hstruct)
  
  if(!disconnected)
  {	
    # print(hypo_mat)  
    
    # print(which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 ))
    
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
    
    # # 	print(hypo_mat)
  }
  
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

  e = edges(graph)
  edge_names = edgeNames(graph)
  nAttrs = list()
  eAttrs = list()
  nAttrs$label = z
  eAttrs$label = e
  names(eAttrs$label) = edge_names
  print(eAttrs)
  
  # set a default color
  nAttrs$fillcolor =  nAttrs$label
  nAttrs$fillcolor[] = 'White'
  
  # use colors defined in tronco$types
  w = unlist(lapply(names(nAttrs$fillcolor), function(x){
    if (x %in% rownames(data$annotations))
      data$types[data$annotations[x,'type'], 'color']
    else
      'White'
    }))
  nAttrs$fillcolor[] = w
  
  # hide node border
  nAttrs$color = nAttrs$fillcolor
  
  # set color for logic nodes
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
  nAttrs$fillcolor[which(w)] = 'orange'
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
  nAttrs$fillcolor[which(w)] = 'green'
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
  nAttrs$fillcolor[which(w)] = 'red'
    
  attrs <- list(node = list(fixedsize = FALSE, fontsize=font, fillcolor='yellow')) 
  
  cur.dev = dev.cur()
  
  pdf(file=paste(name, as.character(disconnected), '.', as.character(pf),'.pdf', sep=''), height=11, width=8.5)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs)
  
  dev.off()
  dev.set(which=cur.dev)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs)
  
}
