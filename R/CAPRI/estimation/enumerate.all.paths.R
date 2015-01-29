#### enumerate.all.paths.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#enumerate all the paths between two nodes of a DAG
#INPUT:
#ancestor.node: first node of the path
#child.node: last node of the path
#parents.pos: topological connections
#RETURN:
#all.paths: vector of all the paths
"enumerate.all.paths" <-
function( ancestor.node, child.node, parents.pos ) {
	#set the initial parents set
	all.paths = parents.pos[[child.node]];
	if(length(all.paths)==1 && all.paths==-1) {
		is.done = TRUE;
	}
	else {
		is.done = FALSE;
	}
	#visit all the nodes in topological order
	while(is.done==FALSE) {
		curr.paths = vector();
		is.done = TRUE;
		for (i in 1:length(all.paths)) {
			curr.new.path = all.paths[[i]];
			curr.new.parents = parents.pos[[curr.new.path[1]]];
			if(length(curr.new.parents)>1 || curr.new.parents!=-1) {
				is.done = FALSE;
				for (j in 1:length(curr.new.parents)) {
					curr.paths[length(curr.paths)+1] = list(c(curr.new.parents[j],curr.new.path));
				}
			}
			else {
				curr.paths[length(curr.paths)+1] = list(curr.new.path);
			}
		}
		all.paths = curr.paths;
	}
	#remove all the paths that are not visiting ancestor.node
	curr.paths = vector();
	for (i in 1:length(all.paths)) {
		curr.result = which(all.paths[[i]]%in%ancestor.node);
		if(length(curr.result)>0) {
			curr.new.path = all.paths[[i]];
			curr.paths[length(curr.paths)+1] = list(c(curr.new.path[which(all.paths[[i]]%in%ancestor.node):length(curr.new.path)],child.node));
		}
	}
	all.paths = unique(curr.paths);
    return(all.paths);
}

#### end of file -- enumerate.all.paths.R
