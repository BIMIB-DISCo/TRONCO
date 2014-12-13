test_types.add <- function(){
	types.add("type", "red")
	checkTrue(any(types[,"color"] == "red"))
	reset()
}

test_types.add_duplicates <- function(){
	types.add("type", "red")
	types.add("type", "red")
	checkTrue(nrow(types) == 1)
	obs <- tryCatch(types.add("type", "red"), warning=conditionMessage)
    checkIdentical("Event type type redefined, now has color: red", obs)
	reset()
}

test_types.add_diff_color <- function(){
	types.add("gain", "red")
	types.add("gain", "blue")
	checkTrue(types[which(types[,"type"] == "gain"), "color"] == "blue")
	reset()
}

test_events.add_duplicates <- function(){
	types.add("gain", "blue")
	events.add("a", "gain", 1)
	events.add("a", "gain", 1)
	checkTrue(nrow(events) == 1)
	reset()
}

test_events.add_same_col <- function(){
	types.add("gain", "blue")
	events.add("a", "gain", 1)
	events.add("b", "gain", 1)
	checkTrue(nrow(events) == 1)
	reset()
}

test_events.add_same_key <- function(){
	types.add("gain", "blue")
	events.add("a", "gain", 1)
	events.add("a", "gain", 2)
	checkTrue(nrow(events) == 1)
	reset()
}

test_events.add <- function(){
	types.add("gain", "blue")
	types.add("loss", "green")	
	events.add("a", "gain", 1)
	events.add("a", "loss", 2)
	checkTrue(nrow(events) == 2)
	reset()
}

