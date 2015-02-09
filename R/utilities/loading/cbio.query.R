"cbio.query" <- function(cbio.study = NA, cbio.dataset = NA, cbio.profile = NA, genes) {

	cat('*** CGDS plugin for cbio query.')
	require("cgdsr")

	if (is.na(cbio.study)) cat('\nAutomatic CBIO study assessment: off')
	else cat(paste('\nAutomatic CBIO study index: ', cbio.study, sep=''))

	if (is.na(cbio.dataset)) cat('\nAutomatic CBIO dataset assessment: off')
	else cat(paste('\nAutomatic CBIO dataset index: ', cbio.dataset, sep=''))

	if (is.na(cbio.profile)) cat('\nAutomatic CBIO profile assessment: off')
	else cat(paste('\nAutomatic CBIO profile index: ', cbio.profile, sep=''))

	mycgds = CGDS("http://www.cbioportal.org/public-portal/")
	
	if (is.na(cbio.study))
	{
		cat("\nAvailable studies at CBIO portal.\n")
		print(getCancerStudies(mycgds)[c("cancer_study_id", "name")])	
		cbio.study <- readline(prompt = "Enter CBIO study id: ")
	}
	
	# Get available case lists (collection of samples) for a given cancer study
	mycancerstudy <- getCancerStudies(mycgds)[cbio.study, ]

	if (is.na(mycancerstudy[1, 1])) 
		stop("Error, CBIO study id invalid. Aborting.")
	
	study.name <- mycancerstudy[1, 1]
	study.ref <- mycancerstudy[1, 2]
	study.syn <- mycancerstudy[1, 3]

	cat(paste("\nCancer codename: ", study.name, sep = ""))
	cat(paste("\nCancer Ref.: ", study.ref, sep = ""))
	cat(paste("\nCancer Syn.: ", study.syn, sep = ""))


	# Get dataset for the study
	if (is.na(cbio.dataset))
	{
		cat("\nAvailable datasets for study.\n")
		print(getCaseLists(mycgds, study.name)[c("case_list_id", "case_list_name")])
		cbio.dataset <- readline(prompt = "Enter study dataset id: ")
	}

	mycaselist = getCaseLists(mycgds, study.name)[cbio.dataset, ]

	# print(mycaselist)
	if (is.na(mycaselist[1, 1])) 
		stop("Error no data for study")

	data.name <- mycaselist[1, 1]
	data.ref <- mycaselist[1, 2]
	data.syn <- mycaselist[1, 3]
	data.id <- mycaselist[1, 4]

	cat(paste("\nData codename: ", data.name, sep = ""))
	cat(paste("\nData Ref.: ", data.ref, sep = ""))
	cat(paste("\nData Syn.: ", data.syn, sep = ""))
	cat(paste("\nData Id.: ", data.id, sep = ""))

	# Get available genetic profiles	
	if (is.na(cbio.profile))
	{
		cat("\nAvailable genetic profiles for selected datasets.\n")
		print(getGeneticProfiles(mycgds, study.name)[c("genetic_profile_id", "genetic_profile_name", "genetic_alteration_type")])
		cbio.profile <- readline(prompt = "Enter genetic profile id: ")
	}
	
	mygeneticprofile = getGeneticProfiles(mycgds, study.name)[cbio.profile, ]

	if (is.na(mycaselist[1, 1])) 
		stop("error no samples for this case")

	samples.name <- mygeneticprofile[1, 1]
	samples.ref <- mygeneticprofile[1, 2]
	samples.syn <- mygeneticprofile[1, 3]
	samples.id <- mygeneticprofile[1, 4]

	cat(paste("\nSamples codename: ", samples.name, sep = ""))
	cat(paste("\nData Ref.: ", samples.ref, sep = ""))
	cat(paste("\nData Syn.: ", samples.syn, sep = ""))
	cat(paste("\nData Id.: ", samples.id, sep = ""))

	cat("\nQuerying the following list of genes: ")
	cat(paste(genes, collapse=','))

	# Get data slices for a specified list of genes, genetic profile and case list
	data <- getProfileData(mycgds, genes, samples.name, data.name)

	# Get clinical data for the case list
	#myclinicaldata = getClinicalData(mycgds,mycaselist)

	# Export
	cat(paste("\nData retrieved: ", nrow(data), " samples, ", ncol(data), " events.", sep = ""))

	ofile <- paste(study.name, data.name, samples.name, "txt", sep = ".")
	write.table(data, file = ofile)
	cat(paste("\nData exported to file: ", ofile, sep = ""))
	
	return(data)
}
