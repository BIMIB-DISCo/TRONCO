test_no_loaded <- function(){
	
	reset()
	obs <- tryCatch(events.add("A", "loss", 1), error=conditionMessage)
    checkIdentical("types  variable not defined!", obs)
    
}
test_shrinkage_value <- function(){
	
	reset()
	data(types)
	data(events)
	data(ov.cgh)
	data.load(ov.cgh)
	obs <- tryCatch(tronco.caprese(data.values, lambda = 2), error=conditionMessage)
    checkIdentical("Lambda coefficient must be in [0:1]!", obs)
    
}