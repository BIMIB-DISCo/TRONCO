#### is.wholenumber.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Check if the given number is a numeric value without decimals
"is.wholenumber" <-
function(x,tol=.Machine$double.eps^0.5) {
	return(abs(x-round(x))<tol);
}

#### end of file -- is.wholenumber.R
