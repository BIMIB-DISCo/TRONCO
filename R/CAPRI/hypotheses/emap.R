#### emap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the position of the event of a given label.event
emap = function( label.event, dataset, annotations ) {
	col.num = -1;
	events.name = "";
	if(!is.null(dataset) && !is.null(annotations)) {
		if(label.event[2]!="*") {
			curr.events = which(annotations[,"event"]==label.event[1] & annotations[,"type"]==label.event[2]);
		}
		else {
			curr.events = which(annotations[,"event"]==label.event[1]);
		}
		if(length(curr.events)>0) {
			events.name = names(curr.events);
			col.num = which(colnames(dataset)%in%events.name);
		}
	}
	else {
		stop("The dataset must be loaded before the hypotheses' definition!",call.=FALSE);
	}
	results = list(col.num=col.num,events.name=events.name);
	return(results);
}

#### end of file -- emap.R
