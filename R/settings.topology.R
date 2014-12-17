#### settings.topology.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Define all the pieces of information needed for the visualization of the reconstructed topology
"settings.topology" <-
function(merged.columns, deleted.column, verbose = FALSE) {
	#define colors and labels to be assigned to each event
	settings$visualization$labels = c();
	settings$visualization$visualized.labels = c();
	settings$visualization$colors = c();
	for(i in 1:nrow(settings$events)) {
		curr.label = paste(toString(search.event(i)$event),toString(search.event(i)$type),sep=":");
		curr.visualized.label = toString(search.event(i)$event);
		curr.color = search.type.info(search.event(i)$type)$color;
		if(is.na(curr.label) || is.na(curr.visualized.label) || is.na(curr.color)) {
			stop(paste("For event \"",search.event(i)$event, "\" no color can be found, please check the types definition!",sep=""),call.=FALSE);
		}
		settings$visualization$labels = c(settings$visualization$labels,curr.label);
		settings$visualization$visualized.labels = c(settings$visualization$visualized.labels,curr.visualized.label);
		settings$visualization$colors = c(settings$visualization$colors,curr.color);
	}
	#events are merged as specified in merged.columns
  	if(!is.na(merged.columns) && length(merged.columns)>0) {
  		for(i in 1:nrow(merged.columns)) {
  			#get the pair of columns to be merged
  			merged.first = merged.columns[i,1];
  			merged.second = merged.columns[i,2];
  			#in labels are saved the labels in the form event_name:type
  			#save the pair of events to be merged
  			settings$visualization$labels[merged.first] = paste(settings$visualization$labels[merged.first],settings$visualization$labels[merged.second],sep = " - ");
  			#in visualized.labels are saved the labels to be assigned to each node in the final plot
  			settings$visualization$visualized.labels[merged.first] = paste(settings$visualization$visualized.labels[merged.first],settings$visualization$visualized.labels[merged.second],sep = " - ");
  			if(verbose) {
  				cat(paste("Column ",toString(merged.first)," of event ",search.event(merged.first)$event," will be merged with column ",toString(merged.second)," of event: ",search.event(merged.second)$event,".", sep = ""));
  				cat("Color lightyellow will be assigned to the merged event.");
			}
			settings$visualization$colors[merged.first] = "lightyellow";
		}
		#if the type merged is not in the list, add it
		if(!search.type("merged")) {
      		settings$types = rbind(settings$types,settings$data.frame(type="merged",color="lightyellow",stringsAsFactors=FALSE));
		}
	}
	#now I can delete any columns present in the vector deleted.column
	if(length(deleted.column) > 0 && !is.na(deleted.column[1])) {
		settings$visualization$labels = settings$visualization$labels[-1*deleted.column];
		settings$visualization$visualized.labels = settings$visualization$visualized.labels[-1*deleted.column];
		settings$visualization$colors = settings$visualization$colors[-1*deleted.column];
		if(verbose) {
      		for(i in 1:length(deleted.column)) {
				cat(paste("Column ", toString(deleted.column[i]), " refferring to event: ",search.event(deleted.column[i])$event, " of type: ",search.event(deleted.column[i])$type," has been removed.",sep=""));
            }
  		}
  	}
  	return(settings);
}

#### end of file -- settings.topology.R
