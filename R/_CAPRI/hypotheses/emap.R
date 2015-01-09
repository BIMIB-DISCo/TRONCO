#### emap.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Return the position of the event of a given label.effect
"emap" <-
function( label.effect, dataset = NA ) {
	if(!is.na(dataset)) {
		col.num = which(names(dataset)==label.effect);
		if(length(col.num)==0) {
			col.num = -1;
		}
	}
	else {
		stop("The dataset must be loaded before the hypotheses' definition!",call.=FALSE);
	}
	return(col.num);
}

#### end of file -- emap.R
