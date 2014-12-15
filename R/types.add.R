#### types.add.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export types.add
#' @title add the name and the color for a new type of mutation (e.g., copy number gain)
#'
#' @description
#' \code{types.add} sets to the global data frame 'settings' the parameter 'types' where the settings for the types of mutations are defined. The types can be added or modified in any order.
#'
#' @details
#' \code{types.add} defines the name and the color for a new type of event (i.e., mutation). If the type was previously defined, its definition is updated to its last definition. A consistency check is performed to ensure that the new type is valid. The types must be defined before the events are loaded.
#' 
#' @param type.name The type name. All type names are strings.
#' @param color.name The type color. All R's standard colors are allowed.
#' 
#' @examples
#' types.add("gain", "red")
#' 
"types.add" <-
function(type.name, color.name) {
	#check if all the input parameters are given
    if(missing(type.name) || missing(color.name)) {
		stop("Missing one or more parameters for the function types.add: types.add(type.name, color.name).",call.=FALSE);
    }
    #otherwise, I can add the new type
    else {
		#if a global types variable is found, the new definition is added in it
		if(exists("settings") && (length(settings$types)>0)) {
			new.type <- settings$types;
			new.type <- rbind(new.type, data.frame(type=type.name,color=color.name,stringsAsFactors=FALSE));
		}
		#otherwise, I create the new definition from scratch
		else {
      		if(!exists("settings")) {
				#create the variable to save the settings
      			settings = list();
      		}
      		new.type <- data.frame(type=type.name,color=color.name,stringsAsFactors=FALSE);
      	}
      	#the user is free to leave spaces between each element in the definition, the definitions file looks more clear this way
      	new.type <- trim.types(new.type);
      	#the check function performs checks for consistency and correctness
		new.type <- check.types(new.type);
		#add the new type to global environment
		cat(paste("Setting color \"",color.name ,"\" for the events of type \"",type.name,"\".\n", sep =""));
		settings$types = new.type;
        assign("settings",settings,envir=.GlobalEnv);
	}
}

#### end of file -- types.add.R
