# commenti

hypothesis.add.group = function(x, FUN, group, dim.min = 2, dim.max = length(group), min.prob = 0, ...) {
	op = deparse(substitute(FUN))

	#print(length(unlist(group)))
	#print(unlist(group))

	effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste(effect, collapse = ", ")

	
	cat("*** Adding Group Hypotheses\n")
	cat('Group:', paste(group, collapse = ", ", sep = ""))
	cat(' Function:', op)
	cat(' Effect:', effect, '\n')
	flush.console()
	
	# group %in% as.events(x)[, 'event']


	if(min.prob > 0)
	{
		cat('\nFiltering genes within the group with alteration frequency below', min.prob, '\n')
		
		temp = events.selection(x, filter.in.names = group)
		
		#show(temp)
		temp = as.alterations(temp)
		temp = events.selection(temp, filter.freq = min.prob)
		
		group = as.genes(temp)
		cat('New group:', paste(group, collapse = ", ", sep = ""), '\n')
	}
	
	ngroup = length(group)
	if (ngroup < 2) {
		warning("No hypothesis can be created for groups with less than 2 elements.")
		return(x)
	}
	
	hom.group = lapply(group, function(g, x) {
		if (nevents(x, genes = g) > 1) 
			T
		else F
	}, x)
	hom.group = group[unlist(hom.group)]

	gene.hom = function(g, h) {
		if (g %in% h) 
		{
			if( any(rowSums(as.gene(x, genes = g)) > 1) ) return(paste0("OR('", g, "')"))
			else return(paste0("XOR('", g, "')"))
		}
		return(paste0("'", g, "'"))
	}

	max.groupsize = min(dim.max, ngroup)
	min.groupsize = max(2, dim.min)
	if(dim.min > dim.max) stop('ERROR - dim.min > dim.max')
	if(min.groupsize > max.groupsize) stop('ERROR - min.groupsize > max.groupsize')
	
	if (length(hom.group) > 0) 
		cat("Genes with multiple events: ", paste(unlist(hom.group), collapse=', ', sep=''), "\n")
	
	error.summary = data.frame()

	# Get an analytical formula... !
	tot.patterns = 0
	for (i in min.groupsize:max.groupsize) tot.patterns = tot.patterns + ncol(combn(unlist(group), i))
	
	# create a progress bar
	cat('Generating ', tot.patterns ,'patterns [size: min =', max.groupsize,' -  max =', max.groupsize, '].\n')
		
	# pb <- txtProgressBar(0, tot.patterns, style = 3)
	flush.console()

	pbPos = 0
	for (i in min.groupsize:max.groupsize) {
		gr = combn(unlist(group), i)
	
		# print(gr)
		
		
		for (j in 1:ncol(gr)) {
			genes = as.list(gr[, j])

			#start the progress bar
			pbPos = pbPos + 1
			# setTxtProgressBar(pb, pbPos)

			hypo.name = paste(unlist(genes), sep = "_", collapse = "_")
			hypo.genes = paste(lapply(genes, function(g, hom.group) {
				gene.hom(g, hom.group)
			}, hom.group), collapse = ", ")

			# print(hypo.genes)
			# print(hom.group)

				hypo.add = paste0("hypothesis.add(x, label.formula = '", op, "_", hypo.name, "', lifted.formula = ", op, "(", hypo.genes, "), ", effect, ")")


				# cat('*** Evaluating ', hypo.add, '\n')
				
				err = tryCatch({
					x = eval(parse(text = hypo.add))
				}, error = function(cond) {
					# print(cond)
					m = paste("Error on", hypo.add, ".\n", cond)
					code = strsplit(as.character(cond), " ")[[1]]
					idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

					return(
						data.frame(
							pattern = paste(unlist(genes), collapse = ", ", sep = ""), 
							error = paste(code[idx.errcode:length(code)], collapse = " ")
							))

				}, warning = function(cond) {
					m = paste("Warning on", hypo.add, ".\n", cond)
					return(genes)
				})
				# Dummy errors detection
				if (!("genotypes" %in% names(err))) 
					error.summary = rbind(error.summary, err)
			}
	}
	# close progress bar
	# close(pb)


	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " genes pattern could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible patterns.\n")

	return(x)
}


hypothesis.add.homologous = function(x, ..., genes = as.genes(x), FUN = "OR") {
	# in questa funzione, per ogni gene che ha più di un tipo di alterazione
	# aggiungo l'OR

	effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste(effect, collapse = ", ")

	hom.group = lapply(genes, function(g, x) {
		if (nevents(x, genes = g) > 1) 
			T
		else F
	}, x)
	hom.group = genes[unlist(hom.group)]

	cat("*** Adding hyoptheses for Homolgous Patterns\n")
	cat('Genes:', paste(hom.group, collapse = ", ", sep = ""))
	cat(' Function:', FUN)
	cat(' Effect:', effect, '\n')
	flush.console()
	
	# create a progress bar
	pb <- txtProgressBar(0, length(hom.group), style = 3)

	error.summary = data.frame()

	for (i in 1:length(hom.group)) {

		#start the progress bar
		setTxtProgressBar(pb, i)

		# Check if the joint probability of homologous events is > 0, if
		# yes the event will be added as 'OR', otherwise 'XOR'
		if( any(
				rowSums(as.gene(x, genes = hom.group[[i]])) > 1) 
			)
		FUN = 'OR'
		else FUN = 'XOR'				

		hypo.add = paste0("hypothesis.add(x, label.formula = '", FUN, "_", hom.group[[i]], "', lifted.formula = ", FUN, "('", hom.group[[i]], "'), ", effect, 
			")")

		err = tryCatch({
			x = eval(parse(text = hypo.add))
		}, error = function(cond) {
			m = paste("Error on", hypo.add, ".\n", cond)
			code = strsplit(as.character(cond), " ")[[1]]
			idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

			return(data.frame(pattern = paste(unlist(hom.group[[i]]), collapse = ", ", sep = ""), error = paste(code[idx.errcode:length(code)], collapse = " ")))

		}, warning = function(cond) {
			m = paste("Warning on", hypo.add, ".\n", cond)
			return(genes)
		})
		# Dummy errors detection
		if (!("genotypes" %in% names(err))) 
			error.summary = rbind(error.summary, err)
	}

	# close progress bar
	close(pb)

	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " patterns could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible gene patterns.\n")


	return(x)

}
