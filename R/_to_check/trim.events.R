#### trim.events.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Trim the spaces at the beginning and in the end of the names of the events variable
"trim.events" <-
function(events) {
	for(i in 1:nrow(events)) {
		for(j in 1:(ncol(events)-1)) {
			events[i,j] <- trim(toString(events[i,j]));
		}
	}
	return(events);
}

#### end of file -- trim.events.R
