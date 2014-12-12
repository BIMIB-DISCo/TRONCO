test_topology <- function(){
	
	reset()
	
	types.add("gain", "red")
	types.add("loss", "blue")
	
	events.add("8q+", "gain", 1)
	events.add("3q+", "gain", 2)
	events.add("5q-", "loss", 3)
	events.add("4q-", "loss", 4)
	events.add("8p-", "loss", 5)
	events.add("1q+", "gain", 6)
	events.add("Xp-", "loss", 7)
	
	data(ov.cgh)
	data.load(ov.cgh)
	
	checkTrue(exists("data.values"))
	checkTrue(ncol(data.values) == 7)
	
	topology <- tronco.caprese(data.values, lambda = 0.5)
	
	adj.matrix <- topology@adj.matrix
	
	# Checks if the expected edges between well known relations are set
	
	checkTrue(adj.matrix[1,2] == 1)
	checkTrue(adj.matrix[3,4] == 1)
	checkTrue(adj.matrix[1,5] == 1)
	checkTrue(adj.matrix[5,7] == 1)
	
	
	set.seed(1234)
	
	topology <- tronco.bootstrap(topology, type = "non-parametric", nboot = 100)
	
	checkEqualsNumeric(topology@confidence$overall.value, 7)
	checkEqualsNumeric(topology@edge.confidence["8q+:gain","3q+:gain"], 0.96)
	checkEqualsNumeric(topology@edge.confidence["5q-:loss","4q-:loss"], 0.45)
	
	
}