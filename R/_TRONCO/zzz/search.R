#### search.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Check if the type with the given tipe.name exists
search.type <- function(type.name) {
	if(exists("settings") && (length(settings$types)>0)) {
		types <- settings$types;
  		for(i in 1:nrow(types)) {
  			if(toString(types[i,"type"])==type.name) {
  				return(TRUE);
  			}
		}
	}
	return(FALSE);
}

# Search for type information for the given type.name
search.type.info <- function(name.type) {
	if(exists("settings") && (length(settings$types)>0)) {
		types <- settings$types;
  		for(i in 1:nrow(types)) {
			if(toString(types[i,"type"])==name.type) {
				return(types[i,]);
			}
		}
    }
    return(NA);
}

# Search for event information for a given column
search.event <- function(column.index) {
	if(exists("settings") && (length(settings$events)>0)) {
		events <- settings$events;
	  	column.index <- as.integer(column.index);
  		for(i in 1:nrow(events)) {
    		if(toString(events[i,"column"])==column.index) {
    			return(events[i,]);
      		}
      	}
    }
  	return(NA);
}

# Search for the event with a given event.name
exists.event <- function(event.name) {
	if(exists("settings") && (length(settings$events)>0)) {
		events <- settings$events;
		for(i in 1:nrow(events)) {
    		if(toString(events[i,"event"])==event.name) {
      			return(events[i,]);
      		}
      	}
    }
  	return(NA);
}

#### end of file -- search.R
