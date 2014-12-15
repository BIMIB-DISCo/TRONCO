#### check.events.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Checks if events variable is correctly defined
"check.events" <-
function(events, types) {
	#looks for duplicates in the declaration
	dup <- c();
	#the flag k is set to TRUE if duplicates by name and type are found
	k <- FALSE;
	#the flag p is set to TRUE if more than one event is assigned to the same column
	p <- FALSE;
	#evaluate the condition of flag k
	i <- 1;
	while(i <= nrow(events)) {
		j <- 1;
		while(j <= (i-1)) {
			#checks if an event is duplicated by name and type
			if((events[j, "event"]==events[i, "event"]) && (events[j, "type"]==events[i, "type"])) {
				dup <- union(dup, c(j));
				k <- TRUE;
			}
			j <- j + 1;
		}
		i <- i + 1;
	}
	if(k) {
		#displays a warning for the duplication
		shown <- c("");
		for(i in 1:length(dup)) {
			if(!any(shown == events[dup[i],"event"])) {
				warning(paste(toString(events[dup[i],"event"])," event already defined!",sep=""),call.=FALSE);
				shown <- union(shown, events[dup[i],"event"]);
			}
		}
	}
	#evaluate the condition of flag p
	dup.col <- c();
	i <- 1;
	while(i <= nrow(events)) {
		j <- 1;
		while(j <= (i-1)) {
			if(events[j,"column"]==events[i,"column"]) {
				dup.col <- union(dup.col,c(j));
				p <- TRUE;
			}
			j <- j + 1;
		}
		i <- i + 1;
	}
	if(p) {
		for(i in 1:length(dup.col)) {
			warning(paste("Event ",toString(events[dup.col[i],"event"])," of type ",toString(events[dup.col[i],"type"])," is in a muitlple defined column and will be discarded", sep=""),call.=FALSE);
		}
 	}
 	#if either p or k are true, remove the duplicated entries
 	if(p || k) {
 		dup <- union(dup,dup.col);
 		events <- events[-1*dup,];
 	}
 	#if the types variable is found (as it should be), check the integrity of each column in the events variable
 	if(exists("settings") && (length(settings$types)>0)) {
 		for(i in 1:nrow(events)) {
 			if(!search.type(toString(events[i,"type"]))) {
 				warning(paste("Type ",toString(events[i,"type"])," can not be found!",sep=""),call.=FALSE);
 			}
 			if(!is.wholenumber(events[i,"column"])) {
 				stop(paste(toString(events[i,"event"])," column must be a valid integer value!",sep=""),call.=FALSE);
 			}
 			if(events[i,"column"]<0) {
 				e <- toString(events[i,"event"]);
 				events <- events[-1*i,];
 				stop(paste(e," column must be an integer positive value: this row will be discarded!",sep=""),call.=FALSE);
      		}
		}
  	}
  	else {
  		stop("The definition of the types can not be found!",call.=FALSE);
  	}
  	#reset the row names of the variable events
  	row.names(events) <- c(1:nrow(events));
  	return(events);
}

#### end of file -- check.events.R
