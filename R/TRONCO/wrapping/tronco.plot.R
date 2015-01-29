#### tronco.plot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


hypotheses.expansion <- function(input_matrix, 
                                 map = list(),
                                 conf_matrix = NULL) {
  
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
      
      # change names in confidence matrix according to hypotesis
      if(!is.null(conf_matrix)) {
        rownames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
        colnames(conf_matrix)[rownames(conf_matrix) == h] = initial_node
      }
      
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
      new_col_name = paste("*", col, sep="_")
      colnames(and_matrix)[pos] = new_col_name
      
      # append a 0 columl to the matrix..
      and_matrix = cbind(and_matrix, matrix(0,nrow = nrow(and_matrix), ncol = 1))
      pos = ncol(and_matrix)
      # and give her the old name
      colnames(and_matrix)[pos] = col
      
      # now do the same to conf_matrix
      if(!is.null(conf_matrix)) {
        conf_matrix = cbind(conf_matrix, conf_matrix[, col])
        pos = ncol(conf_matrix)
        colnames(conf_matrix)[pos] = new_col_name
      }
      
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
    rownames(and_matrix)[pos] = paste0("*_", row)
    and_matrix[paste0("*_", row),row] = 1
  }
  
  # sort col and row (igraph wants the same order)
  and_matrix = and_matrix[,order(colnames(and_matrix))]
  and_matrix = and_matrix[order(rownames(and_matrix)),]
  
  # print(and_matrix)
  if(!is.null(conf_matrix)) {
    return(list(and_matrix, conf_matrix))
  }
  return(and_matrix)
}

###########################
####### TRONCO PLOT #######
###########################


#' @import Rgraphviz
#' @import graph
#' @export tronco.plot
#' @title plot a progression model
#'
#' @description
#' \code{tronco.plot} plots a progression model from a recostructed \code{curr.topology}. 
#' 
#' 
#' @param curr.topology A curr.topology returned by a reconstruction algorithm
#' @param title plot Plot title (default "Progression model x", x reconstruction algorithm)
#' @param title.color color title (default "black")
#' 
#' @param legend bool; show/hide the legend (default is t)
#' @param legend.pos string; legend positioning, available keywords "topleft", "topright","bottomleft" and "bottomright" (default is "bottomright")
#' @param legend.title string; legend title (default is "Legend")
#' 
#' @param legend.columns int; use 1 or 2 columns to plot the legend (default is 1)
#' @param legend.inline bool; print inline legend (default is f)
#' @param legend.coeff double; size of the types label in the legend (default is 1)
#' 
#' @param label.coeff double; size of the events label (default is 1)
#' @param label.color color events label (default "black")
#' @param label.edge.size double; size of the confidence label, when used (default is 12)
#' 
#' @param confidence bool; plot edges according to confidence (default is f)
#' @param node.th.on controls the node thickness, based on the margina probabilty for each event
#' @param node.th; the node thickness, default 2. Ignored if node.th.on is not set.  
#' @examples
#' \dontrun{
#'     types.load("data/types.txt");
#'     events.load("data/events.txt");
#'   	data.load("data/CGH.txt");
#'   	topology <- tronco.caprese();
#'   	tronco.plot(curr.topology, legend.pos = "topleft", legend = TRUE, confidence = TRUE, legend.col = 1, legend.coeff = 0.7, label.edge.size = 10, label.coeff = 0.7);
#' }
tronco.plot = function(x, 
                     fontsize=18, 
                     fontsize.logic=12, 
                     height=1,
                     height.logic=0.5,
                     width=1.5,
                     width.logic=0.5,
                     pf = FALSE, 
                     disconnected=FALSE,
                     fixed.size=FALSE,
                     name=deparse(substitute(capri)),
                     title = paste("Progression model", x$parameters$algorithm), 
                     # title.color = "black", 
                     confidence = FALSE, 
                     legend = TRUE, 
                     legend.title = "Legend", 
                     legend.columns = 1, 
                     legend.inline = FALSE, 
                     legend.pos = "bottomright", 
                     legend.coeff = 1, 
                     # label.coeff = 1, 
                     # label.color = "black", 
                     label.edge.size = 12, 
                     node.th.on = FALSE) {
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
  marginal_p = x$probabilities$probabilities.bic$marginal.probs
  
  
  if (pf) {
    c_matrix = adj.matrix$adj.matrix.pf
    marginal_p = x$probabilities$probabilities.pf$marginal.probs
  }
  
  if (confidence) {
    conf_matrix = if (pf) x$bootstrap$edge.confidence$edge.confidence.pf else x$bootstrap$edge.confidence$edge.confidence.bic
  }
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
  if (!confidence) {
    hypo_mat = hypotheses.expansion(c_matrix, hstruct)
  } else {
    expansion = hypotheses.expansion(c_matrix, hstruct, conf_matrix)
    hypo_mat = expansion[[1]]
    conf_matrix = expansion[[2]]
  }
  
  # remove disconnected nodes
  if(!disconnected) {	
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
  }
  
  attrs = list(node = list(fixedsize = fixed.size))
      
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
  
  # set fontsize
  nAttrs$fontsize = rep(fontsize, length(node_names))
  names(nAttrs$fontsize) = node_names
  
  # set node height
  nAttrs$height = rep(height, length(node_names))
  names(nAttrs$height) = node_names
  
  # set node width
  nAttrs$width = rep(width, length(node_names))
  names(nAttrs$width) = node_names
  
  if (node.th.on) {
    logical_op = list("AND", "OR", "NOT", "XOR", "*")

    # foreach node
    for (node in node_names) {
      prefix = gsub("_.*$", "", node)
      if ( !(prefix %in% logical_op)) {
        increase_coeff = log(marginal_p[node,] + 1.2)
        # print(increase_coeff)
        nAttrs$width[node] = nAttrs$width[node] * increase_coeff
        nAttrs$height[node] = nAttrs$height[node] * increase_coeff
      }
    }
  }
  
  # use colors defined in tronco$types
  w = unlist(lapply(names(nAttrs$fillcolor), function(x){
    if (x %in% rownames(data$annotations))
      data$types[data$annotations[x,'type'], 'color']
    else
      'White'
    }))
  nAttrs$fillcolor[] = w
  
  # hide node border
  nAttrs$color = rep("black", length(node_names))
  names(nAttrs$color) = node_names
  
  # Set shape
  nAttrs$shape = rep("ellipse", length(node_names))
  names(nAttrs$shape) = node_names
  
  # set color, size fo and shape each logic nodes
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
  nAttrs$fillcolor[which(w)] = 'orange'
  nAttrs$label[which(w)] = 'OR'
  nAttrs$shape[which(w)] = 'circle'
  nAttrs$color[which(w)] = 'darkblue'
  nAttrs$fontsize[which(w)] = fontsize.logic
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = width.logic
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
  nAttrs$fillcolor[which(w)] = 'green'
  nAttrs$label[which(w)] = 'AND'
  nAttrs$shape[which(w)] = 'circle'
  nAttrs$color[which(w)] = 'darkblue'
  nAttrs$fontsize[which(w)] = fontsize.logic
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = width.logic
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
  nAttrs$fillcolor[which(w)] = 'red'
  nAttrs$label[which(w)] = 'XOR'
  nAttrs$shape[which(w)] = 'circle'
  nAttrs$color[which(w)] = 'darkblue'
  nAttrs$fontsize[which(w)] = fontsize.logic
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = width.logic
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == '*'
  nAttrs$fillcolor[which(w)] = 'green'
  nAttrs$label[which(w)] = 'AND'
  nAttrs$shape[which(w)] = 'circle'
  nAttrs$color[which(w)] = 'darkblue'
  nAttrs$fontsize[which(w)] = fontsize.logic
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = width.logic
  
  # edges properties
  
  edge_names = edgeNames(graph)
  eAttrs = list()
  
  # set temporary edge shape
  eAttrs$lty = rep("solid", length(edge_names))
  names(eAttrs$lty) = edge_names
  
  #set edge thikness based on prob
  eAttrs$lwd = rep(1, length(edge_names))
  names(eAttrs$lwd) = edge_names
  
  #set edge name based on prob
  eAttrs$label = rep('', length(edge_names))
  names(eAttrs$label) = edge_names
  
  #set arrowdir to forward (default)
  eAttrs$fontsize = rep(label.edge.size, length(edge_names))
  names(eAttrs$fontsize) = edge_names
  
  #set edge color to black (default)
  eAttrs$color = rep('black', length(edge_names))
  names(eAttrs$color) = edge_names
  
  if(pf) {
    # for each edge..
    bic = adj.matrix$adj.matrix.bic
    # print(bic)
    
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      to = edge[2]
      # ..checks if edge is present in BIC
      #print(from)
      #print(to)
      if ( !(from %in% rownames(bic) && to %in% colnames(bic)) ) {
        #print("prima facie!!!")
        #eAttrs$color[e] = 'red'
      }
    }
  }
  
  if(confidence) {
    # for each edge..
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      to = edge[2]
      # ..checks if confidence is available
      if (from %in% rownames(conf_matrix) && to %in% colnames(conf_matrix)) {
        # if confidence > 0..
        # print(paste('from', from, 'to', to, ':', conf_matrix[from, to]))
        if (conf_matrix[from, to] == 1) {
          # ..set edge thickness and label..
          eAttrs$label[e] = '      1'
          eAttrs$lwd[e] = log(150)
        } else if (conf_matrix[from, to] >= 0.01) {
          # ..draw it on the graph..
          eAttrs$label[e] = paste0('      ', substr(conf_matrix[from, to], 2, 4))
          eAttrs$lwd[e] = log(conf_matrix[from, to] * 150)
        } else {
          # ..else set the style of the edge to dashed
          eAttrs$label[e] = "       <.01"
          eAttrs$lwd[e] = log(1.5)
        }
      } else {
        # ..else this edge is located inside to an hypothesis, so no arrow to show
        # eAttrs$dir[e] = 'none'
        eAttrs$color[e] = 'darkblue'
      }
    }
  }
  
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
    if(legend.inline) {
      legend.columns  = length(legend_names);
    }
  }
  
  cur.dev = dev.cur()
  
  pdf(file=paste(name, as.character(disconnected), '.', as.character(pf),'.pdf', sep=''), height=11, width=8.5)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs, main=title)
  # Adds the legend to the plot
  if (legend) {
    legend(legend.pos,
           legend = legend_names,
           title = legend.title,
           bty = 'n',
           cex = legend.coeff,
           pt.cex = 1.5,
           pch = c(19,19),
           ncol = legend.columns,
           col = legend_colors,
           xjust = 1,
           xpd = TRUE)
  }
  
  dev.off()
  dev.set(which=cur.dev)
  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs, main=title)
  # Adds the legend to the plot.
  if (legend) {
    legend(legend.pos,
           legend = legend_names,
           title = legend.title,
           bty = 'n',
           cex = legend.coeff,
           pt.cex = 1.5,
           pch = c(19,19),
           ncol = legend.columns,
           col = legend_colors,
           xjust = 1,
           xpd = TRUE)
  }
  
}
