#### get.lifted.formula.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the adjacency matrix of the formula given the list of edges
"get.lifted.formula" <-
function( lifted.edges ) {
	#structure to save the adjacency matrix
	lifted.adj.matrix = array(0,c(length(unique(c(lifted.edges[,1],lifted.edges[,2]))),length(unique(c(lifted.edges[,1],lifted.edges[,2])))));
	rownames(lifted.adj.matrix) = unique(c(lifted.edges[,1],lifted.edges[,2]));
	colnames(lifted.adj.matrix) = rownames(lifted.adj.matrix);
	#build the matrix given the lifted.edges
	for(i in 1:nrow(lifted.edges)) {
		lifted.adj.matrix[lifted.edges[i,1],lifted.edges[i,2]] = 1;
	}
	return(lifted.adj.matrix);
}

#### end of file -- get.lifted.formula.R
