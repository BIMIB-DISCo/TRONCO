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

###########################
######## HYPO PLOT ########
###########################


hypo.plot = function(x, 
                     font=14, 
                     pf = FALSE, 
                     disconnected=FALSE, 
                     name=deparse(substitute(capri)),
                     # new parameters
                     title = paste("Progression model",curr.topology@algorithm,sep = " "), 
                     title.color = "black", 
                     confidence = FALSE, 
                     legend = TRUE, 
                     legend.title = "Legend", 
                     legend.columns = 1, 
                     legend.inline = FALSE, 
                     legend.pos = "bottomright", 
                     legend.coeff = 1, 
                     label.coeff = 1, 
                     label.color = "black", 
                     label.edge.size = 12, 
                     node.th.on = FALSE, 
                     node.th = 2, 
                     bootstrap="non-parametric") {
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    install.packages('Rgraphviz', dependencies = TRUE)
    library(Rgraphviz)
  }
  
  # Checks if topology exists
  if(missing(x)) {
    stop("Topology missing, usage: hypo.plot(topology, ...", call.=FALSE);
  }
  
  # Want confidence? Boostrap is needed
  if(confidence && !exists('bootstrap', where=x)) {
    stop("To show confidence information, bootstrap execution is needed! See: the function tronco.bootstrap.", call.=FALSE);
  }
  
  # get TRONCO object
  data = x$data
  # print(data)
  
  # get the adjacency matrix
  adj.matrix = x$adj.matrix
  c_matrix = adj.matrix$adj.matrix.bic
  if(pf) c_matrix = adj.matrix$adj.matrix.pf
  # print(c_matrix)
  
  # get algorithm parameters
  parameters = x$parameters
  # print(parameters)
  
  # get hypotheses
  hypotheses = data$hypotheses
  hstruct = NULL
  if (!is.null(hypotheses)) {
    hstruct = hypotheses$hstructure
  }
  
  # expand hypotheses
  hypo_mat = hypotheses.expansion(c_matrix, hstruct)
  
  # remove disconnected nodes
  if(!disconnected) {	
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
  }
  
  hypo_graph = graph.adjacency(hypo_mat)
  v_names = gsub("_.*$", "", V(hypo_graph)$name)
  new_name = list()
  for(v in v_names) {
    if(v %in% rownames(data$annotations)) {
      n = data$annotations[v,"event"]
      new_name = append(new_name, n)
    } else {
      new_name = append(new_name, v)
    }
  }
  V(hypo_graph)$label = new_name
  graph <- igraph.to.graphNEL(hypo_graph)

  node_names = nodes(graph)
  nAttrs = list()
  
  nAttrs$label = V(hypo_graph)$label
  names(nAttrs$label) = node_names
    
  # set a default color
  nAttrs$fillcolor =  rep('White', length(node_names))
  names(nAttrs$fillcolor) = node_names
  
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
  
  # Set shape
  nAttrs$shape = rep("ellipse", length(node_names))
  names(nAttrs$shape) = node_names
  
  # set color, size fo anrd shape each logic nodes
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
  nAttrs$fillcolor[which(w)] = 'orange'
  nAttrs$shape[which(w)] = ''
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
  nAttrs$fillcolor[which(w)] = 'green'
  nAttrs$shape[which(w)] = ''
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
  nAttrs$fillcolor[which(w)] = 'red'
  nAttrs$shape[which(w)] = ''
  
    
  attrs <- list(node = list(fixedsize = FALSE, fontsize=font, fillcolor='yellow')) 
  
  # edges properties
  
  edge_names = edgeNames(graph)
  eAttrs = list()
  
  # set temporary edge name
  #eAttrs$label = rep('    4', length(edge_names))
  #names(eAttrs$label) = edge_names
  
  
  # set temporary edge tick
  #eAttrs$lwd = rep(4, length(edge_names))
  #names(eAttrs$lwd) = edge_names
  
  # set temporary edge shape
  #eAttrs$lty = rep("dashed", length(edge_names))
  #names(eAttrs$lty) = edge_names
  
  # set temporary edge arrow
  eAttrs$arrowhead = rep("open", length(edge_names))
  names(eAttrs$arrowhead) = edge_names
  
  # print(eAttrs)
  
  # create legend
  
  if (legend) {
    legend_colors = data$types[data$types[, "color"] != '#FFFFFF',]
    legend_names = names(legend_colors)
    print(legend_names)
    print(legend_colors)
    if(legend.inline) {
      legend.columns  = length(legend_names);
    }
  }
  
  
  cur.dev = dev.cur()
  
  pdf(file=paste(name, as.character(disconnected), '.', as.character(pf),'.pdf', sep=''), height=11, width=8.5)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs)
  # Adds the legend to the plot
  if (legend) {
    legend(legend.pos,
           legend = legend_names,
           title = legend.title,
           bty = 'n',
           cex = legend.coeff,
           pch = c(19,19),
           #pt.cex = legend.coeff*log(proportion)
           ncol = legend.columns,
           col = legend_colors,
           xjust = 1,
           xpd = TRUE,
           y.intersp = 1.7,
           x.intersp = 1.2)
  }
  
  dev.off()
  dev.set(which=cur.dev)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs)
  # Adds the legend to the plot.
  if (legend) {
    legend(legend.pos,
           legend = legend_names,
           title = legend.title,
           bty = 'n',
           cex = legend.coeff,
           pch = c(19,19),
           #pt.cex = legend.coeff*log(proportion)
           ncol = legend.columns,
           col = legend_colors,
           xjust = 1,
           xpd = TRUE,
           y.intersp = 1.7,
           x.intersp = 1.2)
  }
  
}
