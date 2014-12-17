#### events.load.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export events.load
#' @title load a set of events (e.g., a copy number gain for 8q+) from a file
#'
#' @description
#' \code{events.load} sets to the global data frame 'settings' the parameter 'events' where the settings for the events are defined in a specified file or dataset to be validated.
#'
#' @details
#' \code{events.load} defines a new set of events (i.e., mutations) from a specified file or dataset. The specified file or dataset must be structured as a csv file. All the definitions are a pair of entries specified as follow:
#' 
#' eventName, typeName, columnNumber
#'  ...     , ...     , ...
#' 
#' @param data.input The input file path or a dataset to be validated.
#' 
#' @seealso \code{\link{events.add}}
#' 
"events.load" <-
function(data.input) {
	#check if all the input parameters are given
	if(missing(data.input)) {
		stop("Missing parameter for the function events.load: events.load(data.input).", call.=FALSE);
	}
    #otherwise, I can add the new event
    else {
    	#types must be defined before the events definition
    	if(exists("settings") && (length(settings$types)>0)) {
    		types <- settings$types;
	    	#if the events are given as a data frame
    		if(is.data.frame(data.input)) {
    			events.input <- data.input;
    		}
	    	#if the events are given as a file
    		else {
    			#if the file name is correct, get the file
    			if(file.exists(data.input)) {
    				#set the error message
    				err <- "";
    				message <- "Errors in the definition of the events!";
    				#the definition file may contain errors in the format, hence a try-catch statement is necessary
    				err <- tryCatch(events.input <- suppressWarnings(read.table(data.input,sep = ",",col.names=c("event", "type", "column"),stringsAsFactors=FALSE)),error=function(e) err <- message);
    				if(toString(err)==message) {
    					stop(err,call.=FALSE);
    				}
    			}
	    		#otherwise, show an error message
    			else {
    				stop("The input file can not be found!", call.=FALSE);
    			}
		}
		#if a global events variable is found, the new definition is added in it
    		if(length(settings$events)>0) {
    			#add the new events to the definitions if any
    			if(nrow(events.input)>0) {
    	    		new.event <- rbind(settings$events,events.input);
			}
			else {
				stop("The input file or dataset is empty!",call.=FALSE);
 			}
		}
 		#otherwise, I create the new definition from scratch
		else {
 			#set the new events to the definitions
    			if(nrow(events.input)>0) {
    				new.event <- events.input;
			}
			else {
				stop("The input file or dataset is empty!",call.=FALSE);
    	  		}
    	  	}
    	  	#the user is free to leave spaces between each element in the definition, the definitions file looks more clear this way
    	  	new.event <- trim.events(new.event);
    	  	#the check function performs checks for consistency and correctness
    	  	new.event <- check.events(new.event,types);
		#add the new type to global environment
		for(i in 1:nrow(new.event)) {
			cat(paste("Adding event \"",new.event[i,"event"],"\" of type \"",new.event[i,"type"],"\" (color: \"",toString(search.type.info(new.event[i,"type"])[,"color"]),"\"), dataset column \"",new.event[i,"column"],"\"\n",sep =""));
		}
		settings$events = new.event;
		assign("settings",settings,envir=.GlobalEnv);
    	}
    	else {
    		stop("The types must be defined before the events definition!",call.=FALSE);
  		}
  	}
}

#### end of file -- events.load.R
