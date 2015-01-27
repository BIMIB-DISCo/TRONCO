#### decimal.to.binary.dag.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#convert an integer decimal number to binary
#INPUT:
#num.decimal: decimal integer to be converted
#num.bits: number of bits to be used
#RETURN:
#num.binary: binary conversion of num.decimal
"decimal.to.binary.dag" <-
function(num.decimal, num.bits) {
    #structure where to save the result
    num.binary = rep(0,num.bits);
    #convert the integer decimal number to binary
    pos = 0;
    while(num.decimal>0) {
        #compute the value of the current step
        num.binary[num.bits-pos] = num.decimal %% 2;
        #divide the number by 2 for the next iteration
        num.decimal = num.decimal %/% 2;
        pos = pos + 1;
    }
    return(num.binary);
}

#### end of file -- decimal.to.binary.dag.R
