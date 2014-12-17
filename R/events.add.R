#### events.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export events.add
#' @title add a new event (e.g., a copy number gain for 8q+)
#'
#' @description
#' \code{events.add} sets to the global data frame 'settings' the parameter 'events' where the settings for the events. The events can be added or modified in any order.
#'
#' @details
#' \code{events.add}  defines a new event (i.e., mutation). If the event was previously defined, its definition is updated to its last definition. A consistency check is performed to ensure that the new event is valid. The events must be defined before the events are loaded.
#' 
#' @param event.name The name of the event (e.g., '8q+'). All the names of the events are strings.
#' @param type.name The type name of this event (e.g., 'copy number gain'). The type names must be loaded before adding any event; a consistency check trigger an error if the given type name is not found.
#' @param column.number The column of the dataset to which this event is referring to. The column number must be a valid integer positive value.
#' 
#' @examples
#' types.add("gain", "red")
#' events.add("8q+", "gain", 1)
#' 
"events.add" <-
function(event.name, type.name, column.number = NA) {
	#check if all the input parameters are given
	if(missing(event.name) || missing(type.name) || missing(column.number)) {
		stop("Missing one or more parameters for the function events.add: events.add(event.name, type.name, column.number).",call.=FALSE);
    }
    #otherwise, I can add the new event
    else {
		#types must be defined before the events definition
		if(exists("settings") && (length(settings$types)>0)) {
			types <- settings$types;
			#column number must be a valid integer number
			if(!is.na(column.number) && is.wholenumber(column.number)) {
				#if a global events variable is found, the new definition is added in it
				if(length(settings$events)>0) {
					new.event = settings$events;
					new.event = rbind(new.event,data.frame(event=event.name,type=type.name,column=as.integer(column.number),stringsAsFactors=FALSE));
				}
				else {
					new.event = data.frame(event=event.name,type=type.name,column=as.integer(column.number),stringsAsFactors=FALSE);
				}
				#the user is free to leave spaces between each element in the definition, the definitions file looks more clear this way
      			new.event <- trim.events(new.event);
      			#the check function performs checks for consistency and correctness
				new.event <- check.events(new.event,types);
				#add the new event to global environment
          		cat(paste("Adding event \"",event.name,"\" of type \"",type.name,"\" (color: \"",toString(search.type.info(type.name)[,"color"]),"\"), dataset column \"",column.number,"\"\n", sep =""));
          		settings$events = new.event;
          		assign("settings",settings,envir=.GlobalEnv);
			}
			else {
				stop("Column number must be a valid integer number!",call.=FALSE);
			}
		}
		else {
			stop("The types must be defined before the events definition!",call.=FALSE);
  		}
  	}
}

#### end of file -- events.add.R
