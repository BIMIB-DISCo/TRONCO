#### hypothesis.cycles.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#evaluate cycles involving any hypothesis
#INPUT:
#data: input dataset and its hypotheses
#adj.matrix: adjacency matrix of the topology
#hypotheses.labels: label of the existing hypotheses
"hypothesis.evaluate.cycles" <-
function(data, adj.matrix, hypotheses.labels, weights.matrix) {
	
	#create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
	#create a map structure where to save the atomic events of each hypothesis
	hatomic = new.env(hash=TRUE,parent=emptyenv());
	
	#evaluate all the existing hypotheses
	for (i in 1:length(hypotheses.labels)) {
		
		#evaluate the current hypothesis
		connections = hypothesis.connections(adj.matrix, hypotheses.labels[i]);
		connections = hypothesis.expand.connections(label=hypotheses.labels[i],events=pattern.events(data,hypotheses.labels[i]),incoming=connections$incoming,outgoing=connections$outgoing,hnames=colnames(adj.matrix),hatomic=hatomic,weights.matrix=weights.matrix);
		
		#save the results for the current hypothesis
		ordered.weights = c(ordered.weights,connections$ordered.weights);
		ordered.edges = c(ordered.edges,connections$ordered.edges);
		hatomic = connections$hatomic;
		
	}
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,hatomic=hatomic));
}

#given the adj.matrix, return the incoming and outgoing connections for any hypothesis
#INPUT:
#adj.matrix: adjacency matrix of the topology
#hypotheses.label: label of the hypothesis
"hypothesis.connections" <-
function(adj.matrix, hypotheses.label) {
	incoming = rownames(adj.matrix)[which(adj.matrix[,hypotheses.label]==1)];
	outgoing = colnames(adj.matrix)[which(adj.matrix[hypotheses.label,]==1)];
	connections = list(incoming=incoming,outgoing=outgoing);
	return(connections);
}

#expand and enumerate all the connections incoming or outgoing an hypothesis
#INPUT:
#label: name of the hypothesis
#events: events in the hypothesis
#incoming: incoming connections
#outgoing: outgoing connections
"hypothesis.expand.connections" <-
function(label, events, incoming, outgoing, hnames, hatomic, weights.matrix) {
	
	#create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
    #get the position of the hypothesis
    hypothesis.pos = which(hnames==label);

	#evalutate the incoming and outgoing connections
    curr.edge.pos = 0;
    if(length(incoming)>0) {
		for(i in 1:length(incoming)) {
			ordered.weights = rbind(ordered.weights,weights.matrix[which(hnames==incoming[i]),hypothesis.pos]);
			curr.edge.pos = curr.edge.pos + 1;
        	new.edge <- array(0, c(3,1));
        	new.edge[1,1] = which(hnames==incoming[i]);
        	new.edge[2,1] = hypothesis.pos;
        	new.edge[3,1] = 2;
        	ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
    if(length(outgoing)>0) {
		for(i in 1:length(outgoing)) {
			ordered.weights = rbind(ordered.weights,weights.matrix[hypothesis.pos,which(hnames==outgoing[i])]);
			curr.edge.pos = curr.edge.pos + 1;
        		new.edge <- array(0, c(3,1));
        		new.edge[1,1] = hypothesis.pos;
        		new.edge[2,1] = which(hnames==outgoing[i]);
        		new.edge[3,1] = 1;
        		ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
	
	#add to the map the atomic events of this hypothesis
	hatomic[[toString(hypothesis.pos)]] = which(hnames%in%events);
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,hatomic=hatomic));
}


#### end of file -- hypothesis.cycles.R
