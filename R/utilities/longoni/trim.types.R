#### trim.types.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Trim the spaces at the beginning and in the end of the names of the types variable
"trim.types" <-
function(types) {
	for(i in 1:nrow(types)) {
		for(j in 1:ncol(types)) {
	        types[i,j] <- trim(toString(types[i,j]));
	    }
    }
	return(types);
}

#### end of file -- trim.types.R
