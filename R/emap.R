#### emap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the position of the event of a given label.effect
"emap" <-
function( label.effect ) {
	if(exists("settings") && (length(settings$data.values)>0)) {
		val = which(names(settings$data.values$dataset)==label.effect);
		if(length(val)==0) {
			val = -1;
		}
	}
	else {
		stop("The dataset must be loaded before the hypotheses definition!",call.=FALSE);
	}
	return(val);
}

#### end of file -- emap.R
