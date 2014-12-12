##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2014, Marco Antoniotti, Giulio Caravagna, Alex Graudenzi,        #
# Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis, Giancarlo Mauri, Bud Mishra #
# and Daniele Ramazzotti.                                                        #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the Eclipse Public License v1.0          #
# which accompanies this distribution, and is available at                       #
# http://www.eclipse.org/legal/epl-v10.html and in the include COPYING file      #
#                                                                                #
# Initial contributors:                                                          #
# Giulio Caravagna, Alex Graudenzi, Mattia Longoni and Daniele Ramazzotti.       #
##################################################################################

#' @import Rgraphviz
#' @import graph
#' @export tronco.plot
#' @title plot a progression model
#'
#' @description
#' \code{tronco.plot} plots a progression model from a recostructed \code{topology}. 
#' 
#' 
#' @param topology A topology returned by a reconstruction algorithm
#' @param title plot Plot title (default "Progression model x", x reconstruction algorithm)
#' @param title.color color title (default "black")
#' 
#' @param legend bool; show/hide the legend (default is t)
#' @param legend.pos string; legend positioning, available keywords "topleft", "topright", 
#' "bottom- left" and "bottomright" (default is "bottomright")
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
#' @examples
#' \dontrun{
#' types.load("data/types.txt")
#' events.load("data/events.txt")
#' data.load("data/CGH.txt")
#' topology <- tronco.caprese(data.values)
#' tronco.plot(topology, legend.pos = "topleft", legend = TRUE, confidence = TRUE, 
#' legend.col = 1, legend.coeff = 0.7, label.edge.size = 10, label.coeff = 0.7)
#' }
tronco.plot <- function(topology, title = paste("Progression model", topology@algorithm, sep = " "), 
						title.color = "black",  confidence = FALSE, legend = TRUE, legend.title = "Legend",
            			legend.columns = 1, legend.inline = FALSE, legend.pos = "bottomright", legend.coeff = 1, label.coeff = 1, 
            			label.color = "black", label.edge.size = 12){
  primafacie = FALSE
  ratio = FALSE
  size = FALSE
  if(missing(topology))
    stop("Missing parameter for tronco.plot function: tronco.plot(topology, ...", call. = FALSE)
  if(exists("events") && exists("types") && (length(types) > 0) && (length(events) > 0)){
    
    events <- events
    types <- types
    invalid.events <- invalid.events
    
    # Collects info for the graph building, such as colors and labels the nodes.
    info  <- infoset(invalid.events$merged.events, invalid.events$removed.events)
    
    colors <- info$colors
    n.names <- info$all.vis.labels
    
    # Build the graph object from the adjacency matrix
    if(topology@algorithm == "CAPRESE"){
      adj.matrix <- topology@adj.matrix
      g <- graphAM(adjMat=adj.matrix, edgemode="directed")
    }
    else if(topology@algorithm == "CAPRI"){
      # Creates a graph object using "prima.facie" adj.matrix
      adj.matrix <- topology@adj.matrix
      g <- graphAM(adjMat=adj.matrix, edgemode="directed")
      # Creates a graph object using the "bic" adj.matrix
      adj.matrix <- topology@adj.matrix.bic
      g1 <- graphAM(adjMat=adj.matrix, edgemode="directed")
    }
    
    # If confidece in requested but bootstrap is not performed an error is displayed
    if(confidence && !topology@bootstrap)
      stop("To show confidence information bootstrap execution is needed! see: tronco.bootstrap function!", call. = FALSE)
    
    # Build a list of edges and their thickness
    edw <- c()
    names <- edgeNames(g)
    edge.style <- c()
    if(topology@algorithm == "CAPRESE" && topology@bootstrap){
     	
      low.conficence.edges <- FALSE
      
      for(i in 1:nrow(adj.matrix))
     		for(j in 1:ncol(adj.matrix))
     			if(adj.matrix[i,j] == 1){
     				edw <- c(edw, topology@edge.confidence[i,j])
     			}
      
      # To let tronco.plot draw edges with a grayscale palette call
      # grayscale.color(g,edw) insted of the rep function
      # to assign colors to edge.color variable
      edge.color <- rep("black", length(edgeNames(g)))
      edge.style <-  rep("solid", length(edgeNames(g)))
      ed.name <- edw
      edw.notk <- edw
      
      # Thickness of edges start form 1 but all bootstrap values are
      # lower than one, so a proprional factor is set.
      # An extra 1.5 factor is been required by the thinner edges
      # to be correctly shown
      edw <- edw * 8
      
      low.confidence.edges <- FALSE
      ed.label <- lapply(ed.name, toString)
      for(i in 1:length(ed.label)){
        if(edw.notk[i] == 1)
          ed.label[i] <- substr(ed.label[i],1,4)
        else if(edw.notk[i] == 0){
            ed.label[i] <- ".0"
            edge.style[i]  <- "dashed"
            edw[i] <- 1
            low.confidence.edges <- TRUE
        }
        else
          ed.label[i] <- substr(ed.label[i],2,4)
        
      }
      ed.label <- paste(" ", ed.label, sep="")
      
      names(edw) <- names
      names(ed.name) <- names
      names(ed.label) <- edgeNames(g)
      names(edge.style) <- edgeNames(g)
      
      if(low.confidence.edges)
        cat("Edges with confidence zero will be displayed with dashed edges.\n")
      
    }
    else
      edge.color <- rep("black", length(edgeNames(g)))
    
    # Set a parameter to proportionally set the size of node labels
    max.node.label.len <- max(nchar(colnames(topology@adj.matrix)))
    mean.node.label.len <- mean(nchar(colnames(topology@adj.matrix)))
	  label.node.length <- mean.node.label.len*label.coeff/max.node.label.len
	  
    # Sets parameters for each edge or node
    shape <- rep("ellipse", nrow(adj.matrix))
    arrows <- rep("open", length(edgeNames(g)))
    textColor <- rep(label.color, nrow(adj.matrix))
    label.edge.size <- rep(label.edge.size, length(edgeNames(g)))

                              
    if(any(types$color == label.color))
      warning("Label with same nodes background color!", call. = FALSE)
    
    # Sets names of each pamater
    names(colors) <- colnames(topology@adj.matrix)
    names(shape) <- colnames(topology@adj.matrix)
    names(arrows) <- edgeNames(g)
    names(n.names) <- colnames(topology@adj.matrix)
    names(textColor) <- colnames(topology@adj.matrix)
    names(edge.color) <- edgeNames(g)
    names(label.edge.size) <- edgeNames(g)

    
    if(topology@algorithm == "CAPRESE" && confidence){
    	edgeRenderInfo(g) <- list(lwd = edw, fontsize = label.edge.size, lty = edge.style)
    	eAttrs = list(label = ed.label, arrowhead = arrows, color = edge.color)
    }
    else 
    	eAttrs = list()
    
    if(topology@algorithm == "CAPRI" && primafacie){
      edge.color <- rep("red", length(edgeNames(g)))
      names(edge.color)  <- edgeNames(g)
      for(i in 1:length(edgeNames(g))){
        if(any(edgeNames(g1) == names(edge.color[i])))
          edge.color[i] <- "black"
      }
      eAttrs = list(arrowhead = arrows, color = edge.color)
      
    }else if(topology@algorithm == "CAPRI" && !primafacie)
      g <- g1
    
    
    # Graph building section.
    # Different type of parameters are assigned to the plot in differents ways. 
    # Please check out all the "RGraphviz" manuals on Bioconductor page and 
    # consider the graph drawing documentation on CRAN
    
    nodeRenderInfo(g) <- list(cex = label.node.length, textCol = textColor)
    
    graph.par(list(graph = list(main = title, col.main = title.color)))
    
    nAttrs = list(shape = shape, fillcolor = colors, color = colors, label = n.names)
    
    if(ratio == FALSE && size == FALSE)
      attrs = list()
    else
    	attrs = list(graph = list(ratio = ratio, size = size))
 
    g <- layoutGraph(g, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs)
    
    #Elimino i bordi attorno al plot
    mar <- par("mar")
    par("mar" = c(1, 1, 1, 1) + 0.1)
    
    # Once the graph is ready to be displayed the renderGraph function does the job.
    renderGraph(g)
    
    if(legend){
      
      legend.names <- types$type
      col <- types$color
      
      if(legend.inline)
        legend.columns  <- length(legend.names)
      
      if(legend.columns > length(legend.names))
        stop("Legend must have at the most columns equal to the number of types!")
      
      # Adds the legend to the plot.
      legend(legend.pos,
  	   legend = legend.names,
       title = legend.title,
  	   bty = 'n',
       cex = legend.coeff,
  	   pch = c(19,19),
  	   ncol = legend.columns,
  	   pt.cex = 3*legend.coeff,
  	   col = col,
  	   xjust = 1,
       xpd = TRUE,
  	   y.intersp = 1.7,
  	   x.intersp = 1.2)
      
    }
    
    par("mar" = mar)
    cat("Plot created successfully\n")
  }
  else
    stop("types, events or topology variable not found")
}