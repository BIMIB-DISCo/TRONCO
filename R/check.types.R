#### check.types.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Checks if types variable is correctly defined
"check.types" <-
function(types) {
	#looks for duplicates in the declaration
	dup <- c();
	#the flag k is set to TRUE if duplicates in type names are found
	k <- FALSE;
	#the flag p is set to TRUE if duplicates in type colors are found
	p <- FALSE;
	#verify all the declarations for duplicates type names
	i <- 1;
	while(i <= nrow(types)) {
		j <- 1;
		while(j <= (i-1)) {
			#checks if any type is duplicated by type name
			if((types[j,"type"] == types[i,"type"])) {
				dup <- union(dup, c(j));
				k <- TRUE;
			}
			j <- j + 1;
		}
		i <- i + 1;
  	}
  	#if any duplicate is found, the last one is kept and a warning message is shown
  	if(k) {
  		types.dup <- types;
  		types <- types[-1*dup,];
  		shown <- c("");
  		for(i in 1:length(dup)) {
  			if(!any(shown == types.dup[dup[i],"type"])) {
  				r <- match(types.dup[dup[i],"type"],types[,"type"]);
  				warning(paste("Event type ",toString(types.dup[dup[i], "type"])," redefined with the new color: ",toString(types[r,"color"]),".",sep = ""),call.=FALSE);
  				shown <- union(shown, types.dup[dup[i],"type"]);
  			}
 		}
  	}
	#verify all the declarations for duplicates type colors
  	i <- 1;
  	while(i <= nrow(types)) {
  		j <- 1;
  		while(j <= (i-1)) {
  			#checks if any type is duplicated by type color
  			if((types[j,"color"] == types[i,"color"]) && (types[j,"type"] != types[i,"type"])) {
  				p <- TRUE;
  			}
  			j <- j + 1;
  		}
  		i <- i + 1;
  	}
  	#if any duplicate is found, a warning message is shown
  	if(p) {
  		warning("Multiple events have been associated to the same color.",call.=FALSE);
  	}
  	#reset the row names of the variable types
  	row.names(types) <- c(1:nrow(types));
  	return(types);
}

#### end of file -- check.types.R
