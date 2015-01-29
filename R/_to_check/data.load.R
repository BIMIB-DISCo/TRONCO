#### data.load.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#' @export data.load
#' @title load a dataset, i.e., a binary matrix, from a file or a RData.
#'
#' @description
#' \code{data.load} save the dataset and the number of its atomic events and samples as an input to a global variable 'settings$data.values'.
#'
#' @details
#' \code{data.load} load a dataset and associates its columns to a specified mutational event. Thus, types and events must be define before calling this function in order to perform a consistency check (see \code{types.load}, \code{types.add},  \code{events.load}, \code{events.add} to load/add types/events).
#' 
#' @param data.input The inputs are the file name where the data are saved or a dataset as a RData. 
#' 
"data.load" <-
function(data.input) {
	if(missing(data.input)) {
		stop("Missing parameter for the function data.load: data.load(data.input).",call.=FALSE);
	}
	else if(exists("settings") && length(settings$types)>0 && length(settings$events)>0) {
		#load the dataset either from a file or from a RData
		if(!is.data.frame(data.input)) {
			settings$data.values$dataset = suppressWarnings(read.table(data.input,header=FALSE));
		}
		else {
			settings$data.values$dataset = data.input;
		}
		#the number of elements of settings$events must be enough to make an association to the dataset columns
		if(ncol(settings$data.values$dataset)>nrow(settings$events)) {
			stop("Incomplete events definition, not enough events have been defined!");
		}
		#if an event declaration refers to an out of bound column, an error message is displayed
		for(i in 1:nrow(settings$events)) {
			if(settings$events[i,"column"]>ncol(settings$data.values$dataset)) {
				stop(paste("Event ",toString(events[i,"event"])," is associated to an invalid column number."));
			}
		}
		column.names = c();
		#each event is associated to his column
		for(i in 1:ncol(settings$data.values$dataset)) {
			event = search.event(i);
			column.name = c(paste(toString(event["event"])," (",toString(event["type"]),", column ",toString(event["column"]),")",sep=""));
			column.names = c(column.names,column.name);
		}
		#check the dataset to be correct accordingly the theory of probabilistic causation
		tmp.data = check.dataset(settings$data.values$dataset,FALSE);
		#get all the visualization info to be used by the object of class topology afterwards
		settings = settings.topology(tmp.data$invalid.events$merged.events,tmp.data$invalid.events$removed.events,TRUE);
		settings$data.values$dataset = tmp.data$dataset;
		#set the default values for the hypotheses (not yet defined)
		settings$hypotheses$num.hypotheses = 0;
		if(!is.null(nrow(tmp.data$invalid.events$merged.events))) {
			cat("The available events are:");
			events.label = settings$visualization$visualized.labels;
			out.matrix = settings$events;
			merged.types <- tmpdata$invalid.events$merged.events[,1];
			out.matrix[merged.types,"type"] = "merged";
			out.matrix = out.matrix[-1*tmpdata$invalid.events$removed.events,];
			names = cbind(settings$visualization$visualized.labels);
			out.matrix[,"event"] = names;
			#columns numbers
			out.matrix[, "column"] = cbind(1:nrow(out.matrix));
			row.names(out.matrix) = c(1:nrow(out.matrix));
			print(out.matrix);
		}
		colnames(settings$data.values$dataset) = settings$visualization$labels;
		if(!is.data.frame(data.input)) {
			cat(paste("Data successfully loaded and validated from the dataset saved at ",toString(data.input),".",sep=""));
		}
  		else {
  			cat("Data successfully loaded and validated from the input data.frame.");
  		}
  		settings$data.values$num.events = ncol(settings$data.values$dataset);
  		settings$data.values$num.patients = nrow(settings$data.values$dataset);
  		assign("settings",settings,envir=.GlobalEnv);
  	}
  	else {
  		stop("Either the events or the types have not been defined!");
  	}
}

#### end of file -- data.load.R
