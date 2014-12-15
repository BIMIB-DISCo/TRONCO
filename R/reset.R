#### reset.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export reset
#' @export reset.types
#' @export reset.events
#' @export reset.data.values
#' @export reset.hypotheses
#' @name reset
#' @title reset
#' @description
#' A set of functions to reset types, events, data.values and hypotheses

#' @rdname reset
#' @usage reset()
#' @details \code{reset()} Resets settings$types, settings$events, settings$data.values, settings$visualization and settings$hypotheses
#' @examples
#' reset()
"reset" <-
function() {
	if(exists("settings")) {
		assign("settings",NULL,envir=.GlobalEnv);
	}
}

#' @rdname reset.types
#' @usage reset.types()
#' @details \code{reset.types()} Resets settings$types
#' @examples
#' reset.types()
"reset.types" <-
function() {
	if(exists("settings")) {
		if(length(settings$types)>0) {
			settings$types = NULL;
			assign("settings",settings,envir=.GlobalEnv);
		}
	}
}

#' @rdname reset.events
#' @usage reset.events()
#' @details \code{reset.events()} Resets settings$events
#' @examples
#' reset.events()
"reset.events" <-
function() {
	if(exists("settings")) {
		if(length(settings$events)>0) {
			settings$events = NULL;
			assign("settings",settings,envir=.GlobalEnv);
		}
	}
}

#' @rdname reset.data.values
#' @usage reset.data.values()
#' @details \code{reset.data.values()} Resets settings$data.values
#' @examples
#' reset.data.values()
"reset.data.values" <-
function() {
	if(exists("settings")) {
		if(length(settings$data.values)>0) {
			settings$data.values = NULL;
			assign("settings",settings,envir=.GlobalEnv);
		}
	}
}

#' @rdname reset.visualization
#' @usage reset.visualization()
#' @details \code{reset.visualization()} Resets settings$visualization
#' @examples
#' reset.visualization()
"reset.visualization" <-
function() {
	if(exists("settings")) {
		if(length(settings$visualization)>0) {
			settings$visualization = NULL;
			assign("settings",settings,envir=.GlobalEnv);
		}
	}
}

#' @rdname reset.hypotheses
#' @usage reset.hypotheses()
#' @details \code{reset.hypotheses()} Resets settings$hypotheses
#' @examples
#' reset.hypotheses()
"reset.hypotheses" <-
function() {
	if(exists("settings")) {
		if(length(settings$hypotheses)>0) {
			settings$hypotheses = NULL;
			assign("settings",settings,envir=.GlobalEnv);
		}
	}
}

#### end of file -- reset.R
