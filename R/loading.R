# Transforms a matrix 'geno' in an input list conformant to TRONCO's specifications. If this casting is not possible
# errors are thrown. If no geno.annot/stage.annot parameters are defined, column names of geno are used as event names.
# If these parameters are specified, they should be compliant with TRONCO's input (explain this better)
#
# - geno: a dataframe with the genotypes (this is not constrained to be 0/1)
# - stage.annot: a nrow(geno)x2 dataframe where stage.annot[i,] denotes the stage associated to geno[i,]. 
# - event.type: type label to be used
# - color: R Brewer palette from which random colors are sampled
# Returns: a list of dataframes compliant with TRONCO's specifications.
#' @import RColorBrewer
#' @export
import.genotypes = function(geno, stage.annot = NA, event.type = "variant", color = "Darkgreen") {

	# Avoid malformed datasets
	if (ncol(geno) == 0 || nrow(geno) == 0) 
		stop("Empty genotypes (number of rows/columns 0), will not import.")
	nc = ncol(geno)
	nr = nrow(geno)

	# Gather col/row names
	if (is.null(colnames(geno))) {
		cn = paste0("Gene", 1:ncol(geno))
		warning("Missing column names to identify genes. Will use labels \"Gene1\", \"Gene2\", .....")
	} else cn = colnames(geno)

	# Gather col/row names
	if (is.null(rownames(geno))) {
		rn = paste0("Sample", 1:nrow(geno))
		warning("Missing row names to identify samples. Will use labels \"Sample1\", \"Sample2\",  .....")
	} else rn = rownames(geno)

	x = list()

	# Access keys - G1, G2, ...
	keys = paste0("G", 1:ncol(geno))

	# Genotype matrix
	x$genotypes = as.matrix(geno)
	colnames(x$genotypes) = keys
	rownames(x$genotypes) = rn

	# Create attributes
	x$annotations = matrix(0, nrow = nc, ncol = 2)
	colnames(x$annotations) = c("type", "event")
	rownames(x$annotations) = keys

	x$annotations[, "type"] = event.type
	x$annotations[, "event"] = cn

	# We create a map from types to colors
	x$types = matrix(color, nrow = 1, ncol = 1)
	rownames(x$types) = event.type
	colnames(x$types) = c("color")

	# If the input color is a ColorBrewer scheme	
	#	my.palette = color	
#if(color %in% rownames(brewer.pal.info)) 
#	my.palette = brewer.pal(n=brewer.pal.info[color, 'maxcolors'], name=color); 
#

	is.compliant(x, "import.genotypes: output")

	return(x)
}


#
# import.gistic - Convert a GISTIC score file to TRONCO input.

#' Return the name annotating the dataset, if any.
#'
#' @title as.name
#' @param x A TRONCO compliant dataset.
#' @return The name annotating the dataset, if any.
#' @export as.name


#' @export
import.GISTIC <- function(x, stage.annot = NA) {
  
  if(is.character(x))
  {
    cat('*** Input "x" is a character, interpreting it as a filename to load a table.
Required table format constitent with TCGA data for focal CNAs:
\t- one column for each sample, one row for each gene;
\t- a column Hugo_Symbol with every gene name;
\t- a column Entrez_Gene_Id with every gene\'s Entrez ID.\n')
    
    data = read.table(x, 
                        header = TRUE, 
                        check.names = F,
                        stringsAsFactors = F)
    
    if(any(is.null(colnames(data)))) stop('Input table should have column names.')    
    if(!'Hugo_Symbol' %in% colnames(data)) stop('Missing Hugo_Symbol column!')
    if(!'Entrez_Gene_Id' %in% colnames(data)) stop('Missing Hugo_Symbol column!')
    data$Entrez_Gene_Id = NULL
    rownames(data) = data$Hugo_Symbol
    data$Hugo_Symbol = NULL
    x = t(data)
  }
  
	cat("*** GISTIC input format conversion started.\n")

	# For next operations it is convenient to have everything as 'char' rather than 'int'
	if (typeof(x[, 1]) != typeof("somechar")) {
		cat("Converting input data to character for import speedup.\n")
		rn = rownames(x)
		x = apply(x, 2, as.character)
		rownames(x) = rn
	}

	if (any(is.na(x))) 
		warning("NA entries were replaced with 0s.\n")
	x[is.na(x)] = 0

	cat("Creating ", 4 * (ncol(x)), "events for", ncol(x), "genes \n")

	# gene symbols
	enames <- colnames(x)
	if (is.null(enames)) 
		stop("Error: gistic file has no column names and can not imported, aborting!")

	cat("Extracting \"Homozygous Loss\" events (GISTIC = -2) \n")
	d.homo = x
	d.homo[d.homo != -2] <- 0
	d.homo[d.homo == -2] <- 1

	cat("Extracting \"Heterozygous Loss\" events (GISTIC = -1) \n")
	d.het <- x
	d.het[d.het != -1] <- 0
	d.het[d.het == -1] <- 1

	cat("Extracting \"Low-level Gain\" events (GISTIC = +1) \n")
	d.low <- x
	d.low[d.low != 1] <- 0

	cat("Extracting \"High-level Gain\" events (GISTIC = +2) \n")
	d.high <- x
	d.high[d.high != 2] <- 0
	d.high[d.high == 2] <- 1

	cat("Transforming events in TRONCO data types ..... \n")
	d.homo = trim(import.genotypes(d.homo, event.type = "Homozygous Loss", color = "dodgerblue4"))
	d.het = trim(import.genotypes(d.het, event.type = "Heterozygous Loss", color = "dodgerblue1"))
	d.low = trim(import.genotypes(d.low, event.type = "Low-level Gain", color = "firebrick1"))
	d.high = trim(import.genotypes(d.high, event.type = "High-level Gain", color = "firebrick4"))

	d.cnv.all = ebind(d.homo, d.het, d.low, d.high)

	cat("*** Data extracted, returning only events observed in at least one sample \n", 
      "Number of events: n =", nevents(d.cnv.all), "\n", 
      "Number of genes: |G| =", ngenes(d.cnv.all), "\n",
      "Number of samples: m =", nsamples(d.cnv.all), "\n")


	is.compliant(d.cnv.all, "import.gistic: output")
	return(d.cnv.all)
}


#
# import.mutations - import from MAF file
#
#' @export
"import.MAF" <- function(file, sep = "\t", is.TCGA = TRUE) {
	cat("*** Importing from file: ", file, "\n")
	cat("Loading MAF file ...")
	maf = read.delim(file, comment.char = "#", sep = sep, header = TRUE, stringsAsFactors = FALSE)
	cat("DONE\n")

	#### Auxiliary functions to extract information from the MAF file
	# This is the possibly smallest type of information required to prepare a TRONCO file
# If any necessary information is missing, execution is aborted
variants = function(x) {
		if (!("Variant_Classification" %in% colnames(x))) 
			warning("Missing Variant_Classification flag in MAF file.")

		return(unique(x[, "Variant_Classification"]))
	}

	samples = function(x) {
		if (!("Tumor_Sample_Barcode" %in% colnames(x))) 
			stop("Missing Tumor_Sample_Barcode flag in MAF file - will not import.")

		return(unique(x[, "Tumor_Sample_Barcode"]))
	}

	genes = function(x) {
		if (!("Hugo_Symbol" %in% colnames(x))) 
			stop("Missing Hugo_Symbol flag in MAF file - will not import.")

		return(unique(x[, "Hugo_Symbol"]))
	}

	valid.calls = function(x) {
		if (!("Validation_Status" %in% colnames(x))) 
			warning("Missing Validation_Status flag in MAF file.")
		else return(which(x[, "Validation_Status"] == "Valid"))
	}

	as.TCGA.patients = function(x) {
		samples = samples(x)
		patients = substr(samples, 0, 12)

		return(unique(patients))
	}

	# General report about the mutations stored in this MAF
	cat("*** MAF report: ")
	if (is.TCGA) 
		cat("TCGA=TRUE")

	MAF.variants = variants(maf)
	MAF.samples = samples(maf)
	MAF.genes = genes(maf)

	cat("\nType of annotated mutations: \n")
	print(MAF.variants)

	cat("Number of samples:", length(MAF.samples), "\n")

	# If it is TCGA you should check for multiple samples per patient
	if (is.TCGA) {
		TCGA.patients = as.TCGA.patients(maf)
		cat("[TCGA = TRUE] Number of TCGA patients:", length(TCGA.patients), "\n")

		if (length(TCGA.patients) != length(MAF.samples)) 
			warning("This MAF contains duplicate samples for some patients - use TCGA functions for further information")
	}

	which.valid.calls = valid.calls(maf)
	n.valid = length(which.valid.calls)
	cat("Number of annotated mutations:", nrow(maf), "\n")

	if (("Validation_Status" %in% colnames(maf))) 
		cat("Mutations annotated with \"Valid\" flag (%):", round(n.valid/nrow(maf) * 100, 0), "\n")
	else cat("Mutations annotated with \"Valid\" flag (%): missing flag\n")

	cat("Number of genes (Hugo_Symbol):", length(MAF.genes), "\n")

	cat("Starting conversion from MAF to 0/1 mutation profiles (1 = mutation) :")
	cat(length(MAF.samples), "x", length(MAF.genes), "\n")

	flush.console()
	pb <- txtProgressBar(1, nrow(maf), style = 3)

	# Temporary binary matrix
	binary.mutations = matrix(0, nrow = length(MAF.samples), ncol = length(MAF.genes))

	colnames(binary.mutations) = MAF.genes
	rownames(binary.mutations) = MAF.samples

	for (i in 1:nrow(maf)) {
		setTxtProgressBar(pb, i)
		binary.mutations[maf$Tumor_Sample_Barcode[i], maf$Hugo_Symbol[i]] = 1
	}
	close(pb)

	cat("Starting conversion from MAF to TRONCO data type.\n")
	tronco.data = import.genotypes(binary.mutations, event.type = "Mutation")
	is.compliant(tronco.data)

	return(tronco.data)
}

#' @export
extract.MAF.HuGO.Entrez.map = function(file, sep = "\t") {
	cat("*** Importing from file: ", file, "\n")
	cat("Loading MAF file ...")

	maf = read.delim(file, comment.char = "#", sep = sep, header = TRUE, stringsAsFactors = FALSE)
	cat("DONE\n")

	map = unique(maf[, c("Hugo_Symbol", "Entrez_Gene_Id")])
	map.missing = which(map[, "Entrez_Gene_Id"] == 0)
	map.missing = map[map.missing, "Hugo_Symbol", drop = FALSE]

	if (nrow(map.missing) > 0) 
		warning("The are Hugo_Symbol with Entrez_Gene_Id equal to 0")
	else cat("Map seem consistent (non-zero Entrez_Gene_Id) for", nrow(map), "genes")

	return(map)
}



#' @import cgdsr
#' @export
"cbio.query" <- function(cbio.study = NA, cbio.dataset = NA, cbio.profile = NA, genes) {

	cat("*** CGDS plugin for cbio query.")
	require("cgdsr")

	if (is.na(cbio.study)) 
		cat("\nAutomatic CBIO study assessment: off")
	else cat(paste("\nAutomatic CBIO study index: ", cbio.study, sep = ""))

	if (is.na(cbio.dataset)) 
		cat("\nAutomatic CBIO dataset assessment: off")
	else cat(paste("\nAutomatic CBIO dataset index: ", cbio.dataset, sep = ""))

	if (is.na(cbio.profile)) 
		cat("\nAutomatic CBIO profile assessment: off")
	else cat(paste("\nAutomatic CBIO profile index: ", cbio.profile, sep = ""))

	mycgds = CGDS("http://www.cbioportal.org/public-portal/")

	if (is.na(cbio.study)) {
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
	if (is.na(cbio.dataset)) {
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
	if (is.na(cbio.profile)) {
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
	cat(paste(genes, collapse = ","))

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
