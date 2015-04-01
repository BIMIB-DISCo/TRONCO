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
	matomic = new.env(hash=TRUE,parent=emptyenv());
	
	#create a map structure where to save the hypotheses of each atomic event
	mhypotheses = new.env(hash=TRUE,parent=emptyenv());
	
	#evaluate all the existing hypotheses
	for (i in 1:length(hypotheses.labels)) {
		
		#evaluate the current hypothesis
		connections = hypothesis.connections(adj.matrix, hypotheses.labels[i]);
		connections = hypothesis.expand.connections(label=hypotheses.labels[i],events=pattern.events(data,hypotheses.labels[i]),incoming=connections$incoming,outgoing=connections$outgoing,hnames=colnames(adj.matrix),matomic=matomic,weights.matrix=weights.matrix);
		
		#save the results for the current hypothesis
		ordered.weights = c(ordered.weights,connections$ordered.weights);
		ordered.edges = c(ordered.edges,connections$ordered.edges);
		matomic = connections$matomic;
		
	}
	
	#add to the map the link between atomic to pattern
	for (i in 1:ncol(adj.matrix)) {
		if(!is.null(data$atoms[[colnames(adj.matrix)[i]]])) {
			#add to the map the hypotheses of this atomic event
			mhypotheses[[toString(i)]] = which(colnames(adj.matrix)%in%data$atoms[[colnames(adj.matrix)[i]]]);
		}
	}
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,matomic=matomic,mhypotheses=mhypotheses));
}

#given the adj.matrix, return the incoming and outgoing connections for any hypothesis
#INPUT:
#adj.matrix: adjacency matrix of the topology
#hypotheses.label: label of the hypothesis
"hypothesis.connections" <-
function(adj.matrix, hypotheses.label) {
	# cat('\nhl', hypotheses.label, '\n nomi colonne:\n',
	 # paste(colnames(adj.matrix)), 'asd')
	
	#### EROS FIX
	# print('*** PRE')
	# print(hypotheses.label)
	# print('*** POST')
	# foo = hypotheses.label
	hypotheses.label = hypotheses.label[hypotheses.label %in% rownames(adj.matrix)]
	# print(hypotheses.label)
	# if(length(hypotheses.label) == 0) {
		# print(rownames(adj.matrix))
		# print(foo)
		# print(foo %in% rownames(adj.matrix))
		# }
	
	
	incoming = rownames(adj.matrix)[which(adj.matrix[,hypotheses.label]==1)];
	outgoing = colnames(adj.matrix)[which(adj.matrix[hypotheses.label,]==1)];
	connections = list(incoming=incoming,outgoing=outgoing);
	
	#### EROS FIX
	# print(connections)
	return(connections);
}

#expand and enumerate all the connections incoming or outgoing an hypothesis
#INPUT:
#label: name of the hypothesis
#events: events in the hypothesis
#incoming: incoming connections
#outgoing: outgoing connections
"hypothesis.expand.connections" <-
function(label, events, incoming, outgoing, hnames, matomic, weights.matrix) {
	
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
			new.edge <- array(0, c(2,1));
        		new.edge[1,1] = which(hnames==incoming[i]);
        		new.edge[2,1] = hypothesis.pos;
        		ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
    if(length(outgoing)>0) {
		for(i in 1:length(outgoing)) {
			ordered.weights = rbind(ordered.weights,weights.matrix[hypothesis.pos,which(hnames==outgoing[i])]);
			curr.edge.pos = curr.edge.pos + 1;
        		new.edge <- array(0, c(2,1));
        		new.edge[1,1] = hypothesis.pos;
        		new.edge[2,1] = which(hnames==outgoing[i]);
        		ordered.edges[curr.edge.pos] = list(new.edge);
		}
	}
	
	#### EROS FIX
	# print('****')
	# print(hypothesis.pos)
	# print('hnames')
	# print(hnames)
	# print('events')
	# print(events)
	if(length(hypothesis.pos) > 0)
		matomic[[toString(hypothesis.pos)]] = which(hnames%in%events)	#add to the map the atomic events of this hypothesis
	else
		print('hypothesis.pos == 0!')
	
	#return the results
	return(list(ordered.weights=ordered.weights,ordered.edges=ordered.edges,matomic=matomic));
}

#given the hypotheses and the adj.matrix, return the updated adj.matrix
"hypothesis.adj.matrix" <-
function(hypotheses, adj.matrix) {
	if(!is.na(hypotheses[1])) {
		# set the invalid entries in the adj.matrix
		# hypotheses can not be causing other hypotheses
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix),(ncol(adj.matrix)-hypotheses$num.hypotheses+1):ncol(adj.matrix)] = 0;
		# consider the given hypotheses only toward the specified possible effects
		hypotheses.matrix = array(0,c(hypotheses$num.hypotheses,ncol(adj.matrix)-hypotheses$num.hypotheses));		
		for (i in 1:nrow(hypotheses$hlist)) {
			cause = which(hypotheses$hlist[i,1]==colnames(adj.matrix));
			effect = which(hypotheses$hlist[i,2]==colnames(adj.matrix));
			if(length(cause)>0 && length(effect)>0) {
				hypotheses.matrix[cause-ncol(adj.matrix)+hypotheses$num.hypotheses,effect] = 1;
			}
		}
		adj.matrix[(ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix),1:(ncol(adj.matrix)-hypotheses$num.hypotheses)] = hypotheses.matrix;
		for(j in (ncol(adj.matrix)-hypotheses$num.hypotheses+1):nrow(adj.matrix)) {
			for(k in 1:(ncol(adj.matrix)-hypotheses$num.hypotheses)) {
				if(adj.matrix[j,k] == 0) {
					adj.matrix[k,j] = 0;
				}
			}
		}
	}
	return(adj.matrix);
}


#### end of file -- hypothesis.cycles.R
