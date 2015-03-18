#### tronco.plot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


hypotheses.expansion <- function(input_matrix, 
                                 map = list(),
                                 hidden_and = T,
                                 expand = T,
                                 events = NULL,
                                 conf_matrix = NULL
                                 ) {
  
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  num_hypos = length(map)
  
  # get node list
  node_list <- colnames(input_matrix)
  #print('input matrix')
  #print(input_matrix)
  
  # cut input matrix
  margin = length(node_list) - num_hypos
  hypos_new_name = list()
  
  # da finire!!!
  if(is.vector(events) && F) {
    cat('\n remove event is broken!!!\n')
    print(events)
    min_graph = graph.adjacency(input_matrix)
    graph <- igraph.to.graphNEL(min_graph)
    edge_names = edgeNames(graph)
    print(edge_names)
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      print(edge)
      from = edge[1]
      to = edge[2]
      check_from = any(unlist(strsplit(from, '_')) %in% events)
      check_to = any(unlist(strsplit(to, '_')) %in% events)
      if (!(check_from && check_to)) {
        input_matrix[from, to] = 0
      }
    }
    print(input_matrix)
  }
  
  

  cat('*** Hypos expansion:')
  # check if there are hypotheses
  if (num_hypos == 0 || !expand) {
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
      
      hypos_new_name[initial_node] = h
      
      # recreate lost edge
      for (node in final_node) {
        min_graph <- min_graph + edge(initial_node, node)
      }

      # check if there are edge from atomic to hypo and recreate them
      h_edge_in <- input_matrix[,h]
      in_node <- names(h_edge_in)[which(h_edge_in==1)]
      node_in_hypo = V(hypo_graph)$name
      atomic_node_in_hypo = list()
      for (node in node_in_hypo) {
        if(!is.logic.node(node))
        atomic_node_in_hypo = append(atomic_node_in_hypo, node)
      }

      for (pre in in_node) {
        for (post in atomic_node_in_hypo) {
          min_graph <- min_graph + edge(pre, post)
        }
      }
      
    }
    min_matrix = get.adjacency(min_graph, sparse = F)
    
  }
  
  cat(' done')

  
  
  # now expand the hidden AND
  #print(min_matrix)
  
  if(hidden_and == F) {
    # sort col and row (igraph wants the same order)
    min_matrix = min_matrix[,order(colnames(min_matrix))]
    min_matrix = min_matrix[order(rownames(min_matrix)),]
    
    # print(min_matrix)
    if(!is.null(conf_matrix)) {
      return(list(min_matrix, hypos_new_name, conf_matrix))
    }
    return(list(min_matrix, hypos_new_name))
  }
  
  cat('\n*** Expand hidden and:')
  
  and_matrix = NULL
  to_reconnect = list()
  logical_op = list("AND", "OR", "NOT", "XOR")
#   logical_op = list("", "", "", "")
 
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
  
  cat(' done')
  
  # sort col and row (igraph wants the same order)
  and_matrix = and_matrix[,order(colnames(and_matrix))]
  and_matrix = and_matrix[order(rownames(and_matrix)),]
  
  # print(and_matrix)
  if(!is.null(conf_matrix)) {
    return(list(and_matrix, hypos_new_name, conf_matrix))
  }
  return(list(and_matrix, hypos_new_name))
}

is.logic.node <- function(node) {
  if(substr(node, start=1, stop=3) == 'OR_')
    return(TRUE)
  if(substr(node, start=1, stop=4) == 'XOR_')
    return(TRUE)
  if(substr(node, start=1, stop=4) == 'AND_')
    return(TRUE)
  if(substr(node, start=1, stop=4) == 'NOT_')
    return(TRUE)
  return(FALSE)
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
#'    data.load("data/CGH.txt");
#'    topology <- tronco.caprese();
#'    tronco.plot(curr.topology, legend.pos = "topleft", legend = TRUE, confidence = TRUE, legend.col = 1, legend.coeff = 0.7, label.edge.size = 10, label.coeff = 0.7);
#' }
tronco.plot = function(x, 
                     fontsize=18, 
                     height=1,
                     width=1.5,
                     height.logic = 0.5,
                     pf = FALSE, 
                     disconnected=FALSE,
                     scale.nodes=NA,
                     name=deparse(substitute(capri)),
                     title = paste("Progression model", x$parameters$algorithm),  
                     confidence = FALSE, 
                     legend = TRUE, 
                     legend.cex = 1.0, 
                     label.edge.size = 12, 
                     node.th.on = FALSE, # via
                     hidden.and = T,
                     expand = T,
                     genes = NULL

                     #file = .... # print to pdf  
                     ) 
{
  if (!require(igraph)) {
    install.packages('igraph', dependencies = TRUE)
    library(igraph)
  }
  
  if (!require(Rgraphviz)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
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
  
  logical_op = list("AND", "OR", "NOT", "XOR", "*")
  #logical_op = list("", "", "", "")
 
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
  
  if (all(c_matrix == F)) {
    stop('No edge in adjacency matrix! Nothing to show here.')
  }
  
  if (confidence) {
    conf_matrix = if (pf) x$bootstrap$edge.confidence$edge.confidence.pf else x$bootstrap$edge.confidence$edge.confidence.bic
  }
  #print(c_matrix)
  
  # get algorithm parameters
  parameters = x$parameters
  # print(parameters)
  
  # get hypotheses
  hypotheses = data$hypotheses
  hstruct = NULL
  if (!is.null(hypotheses) && !is.na(hypotheses) ) {
    hstruct = hypotheses$hstructure
  }

  # get event from genes list
  events = NULL
  if (is.vector(genes)) {
    events = unlist(lapply(genes, function(x){names(which(as.events(data)[,'event'] == x))}))
  }
  
  # expand hypotheses
  if (!confidence) {
    expansion = hypotheses.expansion(c_matrix, 
                                     hstruct, 
                                     hidden.and, 
                                     expand, 
                                     events)
    hypo_mat = expansion[[1]]
    hypos_new_name = expansion[[2]]
  } else {
    expansion = hypotheses.expansion(c_matrix, 
                                     hstruct, 
                                     hidden.and, 
                                     expand, 
                                     events, 
                                     conf_matrix)
    hypo_mat = expansion[[1]]
    hypos_new_name = expansion[[2]]
    conf_matrix = expansion[[3]]
  }
  
  # print(hypo_mat)
  
  # remove disconnected nodes
  if(!disconnected) { 
    del = which(rowSums(hypo_mat)+colSums(hypo_mat) == 0 )
    w = !(rownames(hypo_mat) %in% names(del))
    hypo_mat = hypo_mat[w,]
    hypo_mat = hypo_mat[,w]
  }
  
  cat('\n*** Render graphics: ')
  
  attrs = list(node = list())
      
  #print(hypo_mat)
  
  hypo_graph = graph.adjacency(hypo_mat)
  #cat('\n')
  #print(V(hypo_graph)$name[26:30])
  v_names = gsub("_.*$", "", V(hypo_graph)$name)
  if (!expand) {
    v_names = gsub("^[*]_(.+)", "*", V(hypo_graph)$name)
  }
  #print(v_names[26:30])
  new_name = list()
  for(v in v_names) {
    if(v %in% rownames(data$annotations)) {
      n = data$annotations[v,"event"]
      new_name = append(new_name, n)
    } else {
      new_name = append(new_name, v)
    }
  }
    
  #print(V(hypo_graph)$name)
  #print(v_names)
  #print(new_name)
  #print(V(hypo_graph)$label)
  #print(hypo_graph)
  
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
  
  if (!is.na(scale.nodes)) {
    
    # foreach node
    min_p = min(marginal_p)
    max_p = max(marginal_p)
    #print(scale.nodes)

    for (node in node_names) {
      prefix = gsub("_.*$", "", node)
      if ( !(prefix %in% logical_op)) {
        #increase_coeff = sqrt(log(marginal_p[node,] + 1) * scale.nodes)
        #increase_coeff = marginal_p[node,]
        increase_coeff = scale.nodes + (marginal_p[node,] - min_p) / (max_p - min_p)
        #increase_coeff = marginal_p[node,] + ((max_p - marginal_p[node,]) / scale.nodes)
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

  legend_logic = NULL
  
  # set color, size form and shape each logic nodes (if hypos expansion actived)
  node.type = 'box'
  if (expand) {
    
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'OR'
    if (any(w)) {
      legend_logic['Exclusivity (soft)'] = 'orange'
    }
    nAttrs$fillcolor[which(w)] = 'orange'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$color[which(w)] = 'black'
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'AND'
    if (any(w)) {
      legend_logic['Co-occurence'] = 'lightgreen'
    }
    nAttrs$fillcolor[which(w)] = 'lightgreen'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$color[which(w)] = 'black'
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
    
    w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == 'XOR'
    if (any(w)) {
      legend_logic['Exclusivity (hard)'] = 'red'
    }
    nAttrs$fillcolor[which(w)] = 'red'
    nAttrs$label[which(w)] = ''
    nAttrs$shape[which(w)] = node.type
    nAttrs$color[which(w)] = 'black'
    nAttrs$height[which(w)] = height.logic
    nAttrs$width[which(w)] = height.logic
  }
  #print(legend_logic)
  
  w = unlist(nAttrs$label[names(nAttrs$fillcolor)]) == '*'
  if (any(w)) {
      legend_logic['Co-occurence'] = 'lightgreen'
    }
  nAttrs$fillcolor[which(w)] = 'lightgreen'
  nAttrs$label[which(w)] = ''
  nAttrs$shape[which(w)] = node.type
  nAttrs$color[which(w)] = 'black'
  nAttrs$height[which(w)] = height.logic
  nAttrs$width[which(w)] = height.logic
  
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
  
  #set fontsize to label.edge.size (default)
  eAttrs$fontsize = rep(label.edge.size, length(edge_names))
  names(eAttrs$fontsize) = edge_names
  
  #set edge color to black (default)
  eAttrs$color = rep('black', length(edge_names))
  names(eAttrs$color) = edge_names
  
  #record logic edge
  eAttrs$logic = rep(F, length(edge_names))
  names(eAttrs$logic) = edge_names
  
  cat('done')
  
  if(confidence) {

    #print(x$confidence[[3]])
    #print(hypos_new_name)

    # for each edge..
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      to = edge[2]
      # ..checks if confidence is available
      if (from %in% rownames(conf_matrix) && to %in% colnames(conf_matrix)) {
        if(from %in% names(hypos_new_name)){ conf_from = hypos_new_name[[from]] } else { conf_from = from }
        if(to %in% names(hypos_new_name)){ conf_to = hypos_new_name[[to]] } else { conf_to = to }
        # if confidence > 0..
        # print(paste('from', from, 'to', to, ':', conf_matrix[from, to]))
        if (conf_matrix[from, to] == 1) {
          # ..set edge thickness and label..
          eAttrs$label[e] = '                  1'
          eAttrs$lwd[e] = log(150)
        } else if (conf_matrix[from, to] >= 0.01) {
          # ..draw it on the graph..

          eAttrs$label[e] = paste0('                  ', substr(conf_matrix[from, to], 2, 4))
          eAttrs$lwd[e] = log(conf_matrix[from, to] * 150)
        } else {
          # ..else set the style of the edge to dashed
          eAttrs$label[e] = "                   <.01"
          eAttrs$lwd[e] = log(1.5)
        }

        hyper_geom = x$confidence[[3]][conf_from, conf_to]
        if (hyper_geom < 0.01) { hyper_geom = '< .01'} else { hyper_geom = round(hyper_geom, 2)}
        eAttrs$label[e] = paste(eAttrs$label[e], ' / ', hyper_geom)


      } else {
        # ..else this edge is located inside to an hypothesis, so no arrow to show
        eAttrs$logic[e] = T
      }
    }
  }

  # remove arrows from logic node (hidden and)
  for(e in edge_names) {
    edge = unlist(strsplit(e, '~'))
    from = substr(edge[1], start=1, stop=1)
    to = edge[2]
    
    if (from == '*') {
      eAttrs$logic[e] = T
    } 
    
    if (is.logic.node(to)) {
      eAttrs$logic[e] = T
    }
  }
  
  #print(eAttrs$lty)
  
  if(pf) {
    cat('\n*** Add prima facie edges: ')
    # for each edge..
    bic = adj.matrix$adj.matrix.bic
    #print(bic)
    #print('logic edge')
    #print(eAttrs$logic)
    
    #print(rownames(bic))
    for(e in edge_names) {
      edge = unlist(strsplit(e, '~'))
      from = edge[1]
      old_name = hypos_new_name[[from]]
      #cat('\n\nnodo from:', hypos_new_name[[from]])
      if (!is.null(old_name)) {
        from = old_name
      }
      to = edge[2]
      if (substr(to, start=1, stop=1) == '*') {
        to = substr(to, start=3, stop=nchar(to))
      }
      
      # ..checks if edge is present in BIC

      #cat('\nfrom:', from, 'to:', to)
      #cat('\nfrom in bic? ', (from %in% rownames(bic)))
      #cat('\nto in bic? ', (to %in% rownames(bic)))
      #cat('\n!eAttrs$logic[e]', !eAttrs$logic[e])
      # check if edge in BIC (valid only if not logic edge) and 'to' is not a fake and
      if ( (from %in% rownames(bic)) &&
           (to %in% colnames(bic)) &&
           !eAttrs$logic[e] &&
           bic[from, to] == 0
           ) {
        #cat("\nprima facie!!!")
        eAttrs$color[e] = 'red'
      } else {
        #cat('\nno PF!')
      }
    }
    cat('done')
  }
  
  

  plot(graph, nodeAttrs=nAttrs, attrs=attrs, edgeAttrs=eAttrs, main=title)
  
  # Adds the legend to the plot
  if (legend) {
    valid_events = colnames(hypo_mat)[which(colnames(hypo_mat) %in% colnames(c_matrix))]
    legend_names = unique(data$annotations[which(rownames(data$annotations) %in% valid_events), 'type'])
    pt_bg = data$types[legend_names, 'color']
    legend_colors = rep('black', length(legend_names))
    pch = rep(21, length(legend_names))
    
    if (length(legend_logic) > 0) {
      pch = c(pch, 0, 0, rep(22, length(legend_logic)))
      legend_names = c(legend_names, ' ', expression(bold('Patterns')), names(legend_logic))
      legend_colors = c(legend_colors, 'white', 'white', rep('black', length(legend_logic)))
      pt_bg = c(pt_bg, 'white', 'white', legend_logic)  
    }
    
    legend('bottomright',
           legend = legend_names,
           title = expression(bold('Events type')),
           bty = 'n',
           cex = legend.cex,
           pt.cex = 1.5 * legend.cex,
           pch = pch,
           col = legend_colors,
           pt.bg = pt_bg)

    #add thickness legend
    valid_names = node_names
    if(expand) {
      valid_names = node_names[unlist(lapply(node_names, function(x){!is.logic.node(x)}))]
    }
    valid_names = grep('^[*]_(.+)$', valid_names, value = T, invert=T)
    dim = nAttrs$height[valid_names]
    prob = marginal_p[valid_names, ]
    
    min = min(dim)
    p_min = round(min(prob) * 100, 0)
    max = max(dim)
    p_max = round(max(prob) * 100, 0)
    
    
    # This is good only if expand = T   
    # throw away hypotheses - cut marginal_p accordingly
    hypo.names = rownames(as.events(x$data, types='Hypothesis'))
    nonhypo.names = setdiff(rownames(as.events(x$data)), hypo.names)
    
    marginal_p = marginal_p[nonhypo.names, , drop = FALSE]

    # Get label of the (first) event with minimum marginale 
    min.p =   rownames(marginal_p)[which(min(marginal_p) == marginal_p) ]
    label.min = as.events(x$data)[ min.p[1] , , drop = FALSE]

    # Get label of the (first) event with max marginale 
    max.p = rownames(marginal_p)[which(max(marginal_p) == marginal_p) ]
    label.max = as.events(x$data)[ max.p[1] , ,  drop = FALSE]

    # Frequency labels
    min.freq = round(min(marginal_p) * 100, 0)
    max.freq = round(max(marginal_p) * 100, 0)
    
    freq.labels = c( 
      paste0(min.freq, ifelse((min.freq < 10 && max.freq > 9), '%  ', '%'), ' ', label.min[, 'event']),
      paste0(max.freq, '% ', label.max[, 'event'])
    )
  
    stat.pch = c(21, 21)
    pt.bg = c(
        as.colors(x$data)[label.min[, 'type']], 
      as.colors(x$data)[label.max[, 'type']]
      )
    col = c('black', 'black')
        
    # Further stats
    y = delete.type(x$data, 'Hypothesis')
    
    freq.labels = c(freq.labels, 
      ' ',
      expression(bold('Sample size')),
      paste0('n = ', nsamples(y)),
      paste0('m = ', nevents(y)),
      paste0('|G| = ', ngenes(y))     
    ) 
     
    stat.pch = c(stat.pch, 0, 0, 20,  20, 20)
    pt.bg = c(pt.bg, 'white', 'white', rep('black', 3))
    col = c(col, 'white', 'white', rep('black', 3)) 
    
    legend('bottomleft',
           legend = freq.labels,
           title = expression(bold('Statistics')),
           bty = 'n',
           box.lty = 3,
           box.lwd = .3,
           pch = stat.pch,
           pt.cex = 1.5  * legend.cex,
           ncol = 1,
           pt.bg = pt.bg,
           cex = legend.cex,
           col = col)
  }
}
