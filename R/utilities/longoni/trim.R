#### trim.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# Remove any space at the beginning or in the end of a string
"trim" <-
function(string) {
	dim <- as.integer(nchar(string));
	str <- substr(string,1,1);
	#removes any left blank space of the string
	while((str == " ")) {
		string <- substr(string,2,dim);
		dim <- as.integer(nchar(string));
		str <- substr(string,1,1);
  	}
  	dim <- as.integer(nchar(string));
  	str <- substr(string,dim,dim);
  	#removes any right blank space of the string
  	while((str == " ")) {
  		string <- substr(string,1,(dim - 1));
  		dim <- as.integer(nchar(string));
  		str <- substr(string,dim,dim);
  	}
  	return(string);
}

#### end of file -- trim.R
