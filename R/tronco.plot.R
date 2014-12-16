#### tronco.plot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


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
#'   	types.load("data/types.txt");
#'   	events.load("data/events.txt");
#'   	data.load("data/CGH.txt");
#'   	topology <- tronco.caprese();
#'   	tronco.plot(curr.topology, legend.pos = "topleft", legend = TRUE, confidence = TRUE, legend.col = 1, legend.coeff = 0.7, label.edge.size = 10, label.coeff = 0.7);
#' }
tronco.plot <- function(curr.topology = NA, title = paste("Progression model",curr.topology@algorithm,sep = " "), title.color = "black", confidence = FALSE, legend = TRUE, legend.title = "Legend", legend.columns = 1, legend.inline = FALSE, legend.pos = "bottomright", legend.coeff = 1, label.coeff = 1, label.color = "black", label.edge.size = 12, node.th.on = FALSE, node.th = 2, primafacie = FALSE, bootstrap="non-parametric") {
	if(suppressWarnings(is.na(curr.topology))) {
		curr.topology = topology;
	}
	if(missing(curr.topology)) {
  		stop("Missing parameter for the function tronco.plot: tronco.plot(curr.topology, ...", call.=FALSE);
  	}
  	if(confidence && !curr.topology@bootstrap.np && !curr.topology@bootstrap.p) {
  		stop("To show confidence information, bootstrap execution is needed! See: the function tronco.bootstrap.", call.=FALSE);
  	}
  	if(exists("settings") && length(settings$types)>0 && length(settings$events)>0 && length(settings$visualization)>0) {
		types <- settings$types;
  		events <- settings$events;
		n.names <- settings$visualization$visualized.labels;
		colors <- settings$visualization$colors;
		# Build the graph object from the adjacency matrix
		if(curr.topology@algorithm == "CAPRESE") {
			adj.matrix <- curr.topology@adj.matrix;
			g <- graphAM(adjMat=adj.matrix,edgemode="directed");
			edgeNames <- edgeNames(g);
		}
		else if(curr.topology@algorithm == "CAPRI") {
			# Creates a graph object using "prima.facie" adj.matrix
			g <- graphAM(adjMat=curr.topology@pf.adj.matrix,edgemode="directed");
			# Creates a graph object using the "bic" adj.matrix
			g.bic <- graphAM(adjMat=curr.topology@adj.matrix,edgemode="directed");
			adj.matrix <- curr.topology@adj.matrix;
			adj.matrix.primafacie <- curr.topology@pf.adj.matrix;
			edgeNames <- edgeNames(g.bic);
			edgeNames.primafacie <- edgeNames(g);
		}
		# Build a list of edges and their thickness
		edw <- c();
		edge.style <- c();
		if((curr.topology@bootstrap.np || curr.topology@bootstrap.p) && confidence) {
    		low.conficence.edges <- FALSE;
    		if(primafacie) {
    			adj.matrix <- adj.matrix.primafacie;
    			edgeNames <- edgeNames.primafacie;
    			if(bootstrap=="non-parametric" && curr.topology@bootstrap.np) {
    				conf.matrix <- curr.topology@pf.edge.confidence.np;
    			}
    			else if(bootstrap=="parametric" && curr.topology@bootstrap.p) {
    				conf.matrix <- curr.topology@pf.edge.confidence.p;
    			}
    			else {
    				if(bootstrap=="non-parametric") {
    					conf.matrix <- curr.topology@pf.edge.confidence.p;
    				}
    				else if(bootstrap=="parametric") {
    					conf.matrix <- curr.topology@pf.edge.confidence.np;
    				}
    			}
    		}
    		else {
    		if(bootstrap=="non-parametric" && curr.topology@bootstrap.np) {
    			conf.matrix <- curr.topology@edge.confidence.np;
    		}
    		else if(bootstrap=="parametric" && curr.topology@bootstrap.p) {
    			conf.matrix <- curr.topology@edge.confidence.p;
    		}
    		else {
    			if(bootstrap=="non-parametric") {
    				conf.matrix <- curr.topology@edge.confidence.p;
    			}
    			else if(bootstrap=="parametric") {
    				conf.matrix <- curr.topology@edge.confidence.np;
    			}
    		}
    		}
    		for(i in 1:nrow(adj.matrix)) {
    			for(j in 1:ncol(adj.matrix)) {
    				if(adj.matrix[i,j] == 1) {
    					edw <- c(edw, conf.matrix[i,j]);
     				}
     			}
     		}
			# To let tronco.plot draw edges with a grayscale palette call
			# grayscale.color(g,edw) insted of the rep function
			# to assign colors to edge.color variable
			edge.color <- rep("black", length(edgeNames));
			edge.style <- rep("solid", length(edgeNames));
			ed.name <- edw;
			edw.notk <- edw;
			# Thickness of edges start form 1 but all bootstrap values are lower than one, so a proprional factor is set.
			edw <- edw * 8;
			low.confidence.edges <- FALSE;
			ed.label <- lapply(ed.name, toString);
			for(i in 1:length(ed.label)) {
				if(edw.notk[i] == 1) {
					ed.label[i] <- substr(ed.label[i],1,4);
				}
				else if(edw.notk[i] == 0) {
					ed.label[i] <- ".0";
					edge.style[i]  <- "dashed";
					edw[i] <- 1;
					low.confidence.edges <- TRUE;
				}
				else {
					ed.label[i] <- substr(ed.label[i],2,4);
				}
			}
			ed.label <- paste(" ", ed.label, sep="");
			names(edw) <- edgeNames;
			names(ed.name) <- edgeNames;
			names(ed.label) <- edgeNames;
			names(edge.style) <- edgeNames;
			if(low.confidence.edges) {
        		cat("Edges with confidence zero will be displayed as dashed lines.\n");
        	}
		}
    	else {
    		edge.color <- rep("black", length(edgeNames));
    	}
    	# Set a parameter to proportionally set the size of node labels
		max.node.label.len <- max(nchar(colnames(curr.topology@adj.matrix)));
		mean.node.label.len <- mean(nchar(colnames(curr.topology@adj.matrix)));
		label.node.length <- mean.node.label.len*label.coeff/max.node.label.len;
		# Sets parameters for each edge or node
		shape <- rep("ellipse", nrow(adj.matrix));
		arrows <- rep("open", length(edgeNames));
		textColor <- rep(label.color, nrow(adj.matrix));
		label.edge.size <- rep(label.edge.size, length(edgeNames));
		if(any(types$color == label.color)) {
    		warning("Label with same nodes background color!", call.=FALSE);
    	}
    	# Sets names of each parameter
		names(colors) <- colnames(curr.topology@adj.matrix);
		names(shape) <- colnames(curr.topology@adj.matrix);
		names(arrows) <- edgeNames;
		names(n.names) <- colnames(curr.topology@adj.matrix);
		names(textColor) <- colnames(curr.topology@adj.matrix);
		names(edge.color) <- edgeNames;
		names(label.edge.size) <- edgeNames;
		if(node.th.on) {
    		# Sets node thickness
    		node.th <- node.th*curr.topology@marginal.probs;
    		node.height.th <- rep(0.5, length(n.names));
    		node.height.th <- node.height.th*node.th;
    		names(node.height.th) <- colnames(curr.topology@adj.matrix);
    		node.width.th <- rep(0.75, length(n.names));
    		node.width.th <- node.width.th*node.th;
    		names(node.width.th) <- colnames(curr.topology@adj.matrix);
		}
		else {
    		node.height.th <- rep(0.5, length(n.names));
    		node.width.th <- rep(0.75, length(n.names));
    		names(node.height.th) <- colnames(curr.topology@adj.matrix);
    		names(node.width.th) <- colnames(curr.topology@adj.matrix);
		}
		eAttrs = list();
		if(curr.topology@algorithm == "CAPRESE" && confidence) {
    		edgeRenderInfo(g) <- list(lwd = edw, fontsize = label.edge.size, lty = edge.style);
    		eAttrs = list(label = ed.label, arrowhead = arrows, color = edge.color);
		}
		if(curr.topology@algorithm == "CAPRI" && primafacie) {
    		edge.color <- rep("red", length(edgeNames.primafacie));
    		names(edge.color) <- edgeNames.primafacie;
    		for(i in 1:length(edgeNames.primafacie)) {
    			if(any(edgeNames(g.bic) == names(edge.color[i]))) {
    				edge.color[i] <- "black";
    			}
    		}
    		eAttrs = list(arrowhead = arrows, color = edge.color);
    		if(confidence) {
    			edgeRenderInfo(g) <- list(lwd = edw, fontsize = label.edge.size, lty = edge.style);
    			eAttrs = list(label = ed.label, arrowhead = arrows, color = edge.color);
    		}
		}
		else if(curr.topology@algorithm == "CAPRI" && !primafacie) {
    		g <- g.bic;
    		if(confidence) {
    			edgeRenderInfo(g) <- list(lwd = edw, fontsize = label.edge.size, lty = edge.style);
    			eAttrs = list(label = ed.label, arrowhead = arrows, color = edge.color);
    		}
    	}
    	# Graph building section.
		# Different type of parameters are assigned to the plot in differents ways.
		# Please check out all the "RGraphviz" manuals on Bioconductor page and 
		# consider the graph drawing documentation on CRAN
		nodeRenderInfo(g) <- list(cex = label.node.length, textCol = textColor);
		graph.par(list(graph = list(main = title, col.main = title.color)));
		nAttrs = list(shape = shape, fillcolor = colors, color = colors, label = n.names, height = node.height.th, width = node.width.th);
		g <- layoutGraph(g, nodeAttrs = nAttrs, edgeAttrs = eAttrs);
		mar <- par("mar");
		par("mar" = c(1, 1, 1, 1) + 0.1);
		# Once the graph is ready to be displayed the renderGraph function does the job.
		renderGraph(g);
		proportion <- c();
		for(i in 1:nrow(types)) {
    		proportion <- c(proportion, sum(events[,"type"] == toString(types[i,"type"])));
    	}
    	if(legend) {
    		legend.names <- types$type;
    		col <- types$color;
    		if(legend.inline) {
    			legend.columns  <- length(legend.names);
    		}
    		if(legend.columns > length(legend.names)) {
    			stop("Legend must have at the most columns equal to the number of types!");
    		}
    		# Adds the legend to the plot.
    		legend(legend.pos,legend = legend.names,title = legend.title,bty = 'n',cex = legend.coeff,pch = c(19,19),ncol = legend.columns,pt.cex = legend.coeff*log(proportion),col = col,xjust = 1,xpd = TRUE,y.intersp = 1.7,x.intersp = 1.2);
    	}
    	par("mar" = mar);
    	cat("Plot created successfully\n");
	}
	else {
    	stop("Types, events or curr.topology not defined.");
    }
}

#### end of file -- tronco.plot.R
