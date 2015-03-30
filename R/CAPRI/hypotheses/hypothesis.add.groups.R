# commenti

hypothesis.add.group = function(x, FUN, group, dim.max = length(group), min.prob = 0, ...) {
	op = deparse(substitute(FUN))

	#print(length(unlist(group)))
	#print(unlist(group))

	effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
	effect = paste(effect, collapse = ", ")

	cat("*** Adding Group Hypotheses\n")
	cat('Group:', paste(group, collapse = ", ", sep = ""), '\n')
	cat('Function:', op, '\n')
	cat('Effect:', effect, "\n")

	if(min.prob > 0)
	{
		cat('Filtering genes within the group with alteration frequency below', min.prob, '\n')
		
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
			return(paste0("OR('", g, "')"))
		return(paste0("'", g, "'"))
	}


	max.groupsize = min(dim.max, ngroup)
	cat('Maximum pattern size:', max.groupsize, '\n')
	
	if (length(hom.group) > 0) 
		cat("[Functional homologous] Genes with multiple events: ", paste(unlist(hom.group), collaps=', ', sep=''), "\n")
	
	error.summary = data.frame()

	for (i in 2:max.groupsize) {
		gr = combn(unlist(group), i)
		
		# create a progress bar
		cat('Generating hypothesis for', ncol(gr), ' patterns with size', i, '\n')		
		pb <- txtProgressBar(1, ncol(gr) + 1, style = 3)
		pbPos = 2
	
		flush.console()

		for (j in 1:ncol(gr)) {
			genes = as.list(gr[, j])

			#start the progress bar
			setTxtProgressBar(pb, pbPos)
			pbPos = pbPos + 1

			hypo.name = paste(unlist(genes), sep = "_", collapse = "_")
			hypo.genes = paste(lapply(genes, function(g, hom.group) {
				gene.hom(g, hom.group)
			}, hom.group), collapse = ", ")

			

				hypo.add = paste0("hypothesis.add(x, label.formula = '", op, "_", hypo.name, "', lifted.formula = ", op, "(", hypo.genes, "), ", effect, ")")

				# cat('*** Evaluating ', hypo.add, '\n')
				
				err = tryCatch({
					x = eval(parse(text = hypo.add))
				}, error = function(cond) {
					m = paste("Error on", hypo.add, ".\n", cond)
					code = strsplit(as.character(cond), " ")[[1]]
					idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

					return(data.frame(pattern = paste(unlist(genes), collapse = ", ", sep = ""), error = paste(code[idx.errcode:length(code)], collapse = " ")))

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

	}


	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " genes pattern could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible gene patterns.\n")

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

	cat("*** Genes with functional homologous found: ", hom.group, "\n")

	# create a progress bar
	# pb <- txtProgressBar(1, length(hom.group), style = 3)


	error.summary = data.frame()

	for (i in 1:length(hom.group)) {

		#start the progress bar
		# setTxtProgressBar(pb, i)

		hypo.add = paste0("hypothesis.add(x, label.formula = '", FUN, "_", hom.group[[i]], "', lifted.formula = ", FUN, "('", hom.group[[i]], "'), ", effect, 
			")")

		err = tryCatch({
			x = eval(parse(text = hypo.add))
		}, error = function(cond) {
			m = paste("Error on", hypo.add, ".\n", cond)
			code = strsplit(as.character(cond), " ")[[1]]
			idx.errcode = which(code == "[ERR]", arr.ind = TRUE) + 1

			return(data.frame(pattern = paste(unlist(genes), collapse = ", ", sep = ""), error = paste(code[idx.errcode:length(code)], collapse = " ")))

		}, warning = function(cond) {
			m = paste("Warning on", hypo.add, ".\n", cond)
			return(genes)
		})
		# Dummy errors detection
		if (!("genotypes" %in% names(err))) 
			error.summary = rbind(error.summary, err)
	}

	# close progress bar
	# close(pb)

	if (nrow(error.summary) > 0) {
		cat(paste(nrow(error.summary), " patterns could not be added -- showing errors\n", sep = ""))
		print(error.summary)
	} else cat("Hypothesis created for all possible gene patterns.\n")


	return(x)

}
