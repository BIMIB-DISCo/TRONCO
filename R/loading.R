#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Import a matrix of 0/1 alterations as a TRONCO compliant dataset. Input "geno" can be either a dataframe or
#' a file name. In any case the dataframe or the table stored in the file must have a column for each altered
#' gene and a rows for each sample. Colnames will be used to determine gene names, if data is loaded from
#' file the first column will be assigned as rownames.
#' 
#' @title import.genotypes
#' @param geno Either a dataframe or a filename
#' @param event.type Any 1 in "geno" will be interpreted as a an observed alteration labeled with type "event.type"
#' @param color This is the color used for visualization of events labeled as of "event.type"
#' @return A TRONCO compliant dataset
#' @export import.genotypes
#'
import.genotypes <- function(geno, event.type = "variant", color = "Darkgreen") {
    if (!(is.data.frame(geno) || is.matrix(geno)) && is.character(geno)) {
        cat('*** Input "geno" is a character, interpreting it as a filename to load a table.
             Required table format:
             \t- one column for each gene, one row for each gene;
             \t- colnames/rownames properly defined.\n')

        data =
            read.table(geno,
                       header = TRUE,
                       check.names = F,
                       stringsAsFactors = F)

        if (any(is.null(colnames(data))))
            stop('Input table should have column names.')
        rownames(data) = data[, 1]
        data[, 1] = NULL
        geno = data
    }

    ## Avoid malformed datasets

    if (ncol(geno) == 0 || nrow(geno) == 0)
        stop("Empty genotypes (number of rows/columns 0), will not import.")
    nc = ncol(geno)
    nr = nrow(geno)

    ## Gather col/row names

    if (is.null(colnames(geno))) {
        cn = paste0("Gene", 1:ncol(geno))
        warning("Missing column names to identify genes. Will use labels \"Gene1\", \"Gene2\", .....")
    } else
        cn = colnames(geno)

    ## Gather col/row names

    if (is.null(rownames(geno))) {
        rn = paste0("Sample", 1:nrow(geno))
        warning("Missing row names to identify samples. Will use labels \"Sample1\", \"Sample2\",  .....")
    } else rn = rownames(geno)

    x = list()

    ## Access keys - G1, G2, ...

    keys = paste0("G", 1:ncol(geno))

    ## Genotype matrix

    x$genotypes = as.matrix(geno)
    colnames(x$genotypes) = keys
    rownames(x$genotypes) = rn

    ## Create attributes

    x$annotations = matrix(0, nrow = nc, ncol = 2)
    colnames(x$annotations) = c("type", "event")
    rownames(x$annotations) = keys

    x$annotations[, "type"] = event.type
    x$annotations[, "event"] = cn

    ## We create a map from types to colors

    x$types = matrix(color, nrow = 1, ncol = 1)
    rownames(x$types) = event.type
    colnames(x$types) = c("color")

    is.compliant(x, "import.genotypes: output")

    return(x)
}


#' Transform GISTIC scores for CNAs in a TRONCO compliant object. Input can be either a matrix, with columns
#' for each altered gene and rows for each sample; in this case colnames/rownames mut be provided. If input
#' is a character an attempt to load a table from file is performed. In this case the input table format
#' should be constitent with TCGA data for focal CNA; there should hence be: one column for each sample,
#' one row for each gene, a column Hugo_Symbol with every gene name and a column Entrez_Gene_Id with every
#'  gene\'s Entrez ID. A valid GISTIC score should be any value of: "Homozygous Loss" (-2), "Heterozygous
#'  Loss" (-1), "Low-level Gain" (+1), "High-level Gain" (+2).
#'
#' @examples
#' data(gistic)
#' gistic = import.GISTIC(gistic)
#' gistic = annotate.description(gistic, 'Example GISTIC')
#' oncoprint(gistic)
#' 
#' @title import.GISTIC
#' @param x Either a dataframe or a filename
#' @return A TRONCO compliant representation of the input CNAs.
#' @export import.GISTIC
#' 
import.GISTIC <- function(x, filter.genes = NULL, filter.samples = NULL) {

    if (!(is.data.frame(x) || is.matrix(x)) && is.character(x)) {
        cat('*** Input "x" is a character, interpreting it as a filename to load a table.
            Required table format constitent with TCGA data for focal CNAs:
            \t- one column for each sample, one row for each gene;
            \t- a column Hugo_Symbol with every gene name;
            \t- a column Entrez_Gene_Id with every gene\'s Entrez ID.\n')

        data =
            read.table(x,
                       header = TRUE,
                       check.names = F,
                       stringsAsFactors = F)

        cat('Data loaded.\n')

        if (any(is.null(colnames(data)))) {
            stop('Input table should have column names.')
        }
        if (!'Hugo_Symbol' %in% colnames(data)) {
            stop('Missing Hugo_Symbol column!')
        }
        if (!'Entrez_Gene_Id' %in% colnames(data)) {
            stop('Missing Hugo_Symbol column!')
        }
        data$Entrez_Gene_Id = NULL
        rownames(data) = data$Hugo_Symbol
        data$Hugo_Symbol = NULL
        x = t(data)
    }

    if(is.null(filter.genes) && is.null(filter.samples)) cat('*** Using full GISTIC: #dim ', nrow(x), ' x ', ncol(x), '\n' ) 
    else cat('*** Filtering full GISTIC: #dim ', nrow(x), ' x ', ncol(x), '\n' ) 

    if(!is.null(filter.genes))
    {
        if(!is.vector(filter.genes))
            stop('filter.genes - should be vector')

        x = x[, which(colnames(x) %in% filter.genes), drop = F ]
        cat('*** Using reduced GISTIC: #dim ', nrow(x), ' x ', ncol(x), '\n' ) 
  }

    if(!is.null(filter.samples))
    {
       if(!is.vector(filter.samples))
            stop('filter.samples - should be vectors')

         x = x[which(rownames(x) %in% filter.samples), ]
        cat('*** Using reduced GISTIC: #dim ', nrow(x), ' x ', ncol(x), '\n' ) 
    }

        

    cat("*** GISTIC input format conversion started.\n")

    ## For next operations it is convenient to have everything as
    ## 'char' rather than 'int'

    if (typeof(x[, 1]) != typeof("somechar")) {
        cat("Converting input data to character for import speedup.\n")
        rn = rownames(x)
        x = apply(x, 2, as.character)
        rownames(x) = rn
    }

    if (any(is.na(x)))  {
        warning("NA entries were replaced with 0s.\n")
    }
    x[is.na(x)] = 0

    cat("Creating ", 4 * (ncol(x)), "events for", ncol(x), "genes \n")

                                        # gene symbols
    enames <- colnames(x)
    if (is.null(enames)) {
        stop("Error: gistic file has no column names and can not imported, aborting!")
    }

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


#' Import mutation profiles from a Manual Annotation Format (MAF) file. All mutations are aggregated as a
#' unique event type labeled "Mutation" and assigned a color according to the default of function
#' \code{import.genotypes}. If this is a TCGA MAF file check for multiple samples per patient is performed
#' and a warning is raised if these occurr.
#'
#' @examples
#' data(maf)
#' mutations = import.MAF(maf)
#' mutations = annotate.description(mutations, 'Example MAF')
#' mutations = TCGA.shorten.barcodes(mutations)
#' oncoprint(mutations)
#' 
#' @title import.MAF
#' @param file  MAF filename
#' @param sep MAF separator, default \'\\t\'
#' @param is.TCGA TRUE if this MAF is from TCGA; thus its sample codenames can be interpreted
#' @return A TRONCO compliant representation of the input MAF
#' @export import.MAF
#' 
import.MAF <- function(file, sep = '\t', is.TCGA = TRUE, filter.fun = NULL) {

    if (!(is.data.frame(file) || is.matrix(file)) && is.character(file)) {
        cat("*** Importing from file: ", file, "\n")
        cat("Loading MAF file ...")
        maf = read.delim(file, comment.char = "#", sep = sep, header = TRUE, stringsAsFactors = FALSE)
        cat("DONE\n")
    } else {
        cat("*** Importing from dataframe\n")
        cat("Loading MAF dataframe ...")
        maf = file
        cat("DONE\n")
    }

    if(is.null(filter.fun)) cat('*** Using full MAF: #entries ', nrow(maf), '\n')
    else{
    	if(!is.function(filter.fun))
    		stop('filter.fun - should be a function')

		cat('*** Filtering full MAF: #entries ', nrow(maf), '\n' ) 
		maf = maf[apply(maf, 1, filter.fun), ]
		cat('*** Using reduced MAF: #entries ', nrow(maf), '\n' ) 
    }


    ## Auxiliary functions to extract information from the MAF file
    ## This is the possibly smallest type of information required to prepare a TRONCO file
    ## If any necessary information is missing, execution is aborted

    variants <- function(x) {
        if (!("Variant_Classification" %in% colnames(x)))
            warning("Missing Variant_Classification flag in MAF file.")

        return(unique(x[, "Variant_Classification"]))
    }

    samples <- function(x) {
        if (!("Tumor_Sample_Barcode" %in% colnames(x)))
            stop("Missing Tumor_Sample_Barcode flag in MAF file - will not import.")

        return(unique(x[, "Tumor_Sample_Barcode"]))
    }

    genes <- function(x) {
        if (!("Hugo_Symbol" %in% colnames(x)))
            stop("Missing Hugo_Symbol flag in MAF file - will not import.")

        return(unique(x[, "Hugo_Symbol"]))
    }

    valid.calls <- function(x) {
        if (!("Validation_Status" %in% colnames(x)))
            warning("Missing Validation_Status flag in MAF file.")
        else return(which(x[, "Validation_Status"] == "Valid"))
    }

    as.TCGA.patients <- function(x) {
        samples = samples(x)
        patients = substr(samples, 0, 12)

        return(unique(patients))
    }

    ## General report about the mutations stored in this MAF

    cat("*** MAF report: ")
    if (is.TCGA) {
        cat("TCGA=TRUE")
    }

    MAF.variants = variants(maf)
    MAF.samples = samples(maf)
    MAF.genes = genes(maf)

    cat("\nType of annotated mutations: \n")
    print(MAF.variants)

    cat("Number of samples:", length(MAF.samples), "\n")

    ## If it is TCGA you should check for multiple samples per patient

    if (is.TCGA) {
        TCGA.patients = as.TCGA.patients(maf)
        cat("[TCGA = TRUE] Number of TCGA patients:", length(TCGA.patients), "\n")

        if (length(TCGA.patients) != length(MAF.samples)) {
            warning("This MAF contains duplicate samples for some patients - use TCGA functions for further information")
        }
    }

    which.valid.calls = valid.calls(maf)
    n.valid = length(which.valid.calls)
    cat("Number of annotated mutations:", nrow(maf), "\n")

    if (("Validation_Status" %in% colnames(maf))) {
        cat("Mutations annotated with \"Valid\" flag (%):", round(n.valid/nrow(maf) * 100, 0), "\n")
    } else {
        cat("Mutations annotated with \"Valid\" flag (%): missing flag\n")
    }
    cat("Number of genes (Hugo_Symbol):", length(MAF.genes), "\n")
    cat("Starting conversion from MAF to 0/1 mutation profiles (1 = mutation) :")
    cat(length(MAF.samples), "x", length(MAF.genes), "\n")

    flush.console()
    pb <- txtProgressBar(1, nrow(maf), style = 3)

    ## Temporary binary matrix

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


#' Extract a map Hugo_Symbol -> Entrez_Gene_Id from a MAF input file. If some genes map to ID 0
#' a warning is raised.
#' 
#' @title extract.MAF.HuGO.Entrez.map
#' @param file  MAF filename
#' @param sep MAF separator, default \'\\t\'
#' @return A mapHugo_Symbol -> Entrez_Gene_Id.
#' @export extract.MAF.HuGO.Entrez.map
#' 
extract.MAF.HuGO.Entrez.map <- function(file, sep = "\t") {
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


#' Wrapper for the CGDS package to query the Cbio portal. This can work either automatically, if one
#' sets \code{cbio.study}, \code{cbio.dataset} or \code{cbio.profile}, or interactively otherwise. A
#' list of genes to query with less than 900 entries should be provided. This function returns a list
#' with two dataframe: the gentic profile required and clinical data for the Cbio study. Output is also
#' saved to disk as Rdata file. See also http://www.cbioportal.org.
#'
#' @title cbio.query
#' @param cbio.study  Cbio study ID
#' @param cbio.dataset Cbio dataset ID
#' @param cbio.profile Cbio genetic profile ID
#' @param genes A list of < 900 genes to query
#' @return A list with two dataframe: the gentic profile required and clinical data for the Cbio study.
#' @export cbio.query
#' @importFrom cgdsr CGDS getCancerStudies getCaseLists getGeneticProfiles getProfileData getClinicalData
#' 
cbio.query <- function(cbio.study = NA, cbio.dataset = NA, cbio.profile = NA, genes) {
    cat("*** CGDS plugin for Cbio query.\n")
    ## require("cgdsr")

    if (is.null(genes) || is.na(genes) || length(genes) == 0) {
        stop('Empty list of genes to query')
    }
    if (length(genes) > 900) {
        stop('URL with more than 900 genes will not be accepted, please split it.')
    }

    if (is.na(cbio.study)) {
        cat("\nAutomatic CBIO study assessment: off")
    }else {
        cat(paste("\nAutomatic CBIO study index: ", cbio.study, sep = ""))
    }

    if (is.na(cbio.dataset)) {
        cat("\nAutomatic CBIO dataset assessment: off")
    } else {
        cat(paste("\nAutomatic CBIO dataset index: ", cbio.dataset, sep = ""))
    }
    if (is.na(cbio.profile)) {
        cat("\nAutomatic CBIO profile assessment: off")
    } else {
        cat(paste("\nAutomatic CBIO profile index: ", cbio.profile, sep = ""))
    }
    mycgds = CGDS("http://www.cbioportal.org/public-portal/")


    cs = getCancerStudies(mycgds)
    if (is.na(cbio.study)) {
        cat("\nAvailable studies at CBIO portal.\n")
        print(cs[c("cancer_study_id", "name")])

        repeat{
            cbio.study <- readline(prompt = "Enter CBIO study id: ")
            if (cbio.study %in% cs$cancer_study_id)
                break
        }
    }

    ## Get available case lists (collection of samples) for a given
    ## cancer study

    mycancerstudy <- cbio.study

    if (is.na(mycancerstudy)) {
        stop("CBIO study id invalid. Aborting.")
    }

    study <- cs[cs$cancer_study_id == cbio.study, , drop = F]

    cat(paste("\nCancer codename: ", study[, 1], sep = ""))
    cat(paste("\nCancer Ref.: ", study[, 2], sep = ""))
    cat(paste("\nCancer Syn.: ", study[, 3], sep = ""))

    cutdescr = function(x, n)
        {
            x[, ncol(x)] = ifelse(
                 nchar(x[, ncol(x)]) > n,
                 paste0(substr(x[, ncol(x)], 1, n), '....'),
                 x[, ncol(x)])
            return(x)
        }

    ## Get dataset for the study

    csl = getCaseLists(mycgds, cbio.study)
    if (is.na(cbio.dataset)) {
        cat("\nAvailable datasets for study:", cbio.study, "\n")
        print(cutdescr(csl[c("case_list_id",  "case_list_description")], 90))

        repeat{
            cbio.dataset <- readline(prompt = "Enter study dataset id: ")
            if (cbio.dataset %in% csl$case_list_id) break
        }
    }

    caselist =  csl[csl$case_list_id == cbio.dataset, , drop = F]

    if (any(is.na(caselist))) stop("No data for selected study. Aborting.")


    cat(paste("\nData codename: ", caselist[, 1], sep = ""))
    cat(paste("\nData Ref.: ", caselist[, 2], sep = ""))
    cat(paste("\nData Syn.: ", caselist[, 3], sep = ""))

    ## Get available genetic profiles

    gp = getGeneticProfiles(mycgds, cbio.study)
    if (is.na(cbio.profile)) {
        cat("\nAvailable genetic profiles for selected datasets.\n")
        print(cutdescr(gp[c("genetic_profile_id", "genetic_profile_description")], 90))

        repeat{
            cbio.profile <- readline(prompt = "Enter genetic profile id: ")
            if (cbio.profile %in% gp$genetic_profile_id) break
        }
    }

    profile = gp[gp$genetic_profile_id == cbio.profile, , drop = F]

    if (any(is.na(cbio.profile))) {
        stop("No samples for this profile. Aborting")
    }

    samples.name <- profile[1, 1]
    samples.ref <- profile[1, 2]
    samples.syn <- profile[1, 3]
    samples.id <- profile[1, 4]

    cat(paste("\nSamples codename: ", samples.name, sep = ""))
    cat(paste("\nData Ref.: ", samples.ref, sep = ""))
    cat(paste("\nData Syn.: ", samples.syn, sep = ""))

    cat("\n\nQuerying the following list of genes: ")
    cat(paste(genes, collapse = ", "), '\n')

    ## Get data slices for a specified list of genes, genetic profile
    ## and case list

    data <- getProfileData(mycgds, genes, samples.name, cbio.dataset)
    rownames(data) = gsub('\\.', '-', rownames(data))
    cat('Symbol \".\" was replaced with "-" in sample IDs.\n')

    ## Export

    cat(paste("\nData retrieved: ", nrow(data), " samples, ", ncol(data), " genes.", sep = ""))

    ## Get clinical data for the case list

    cat("\nRetrieved also clinical data for samples:", cbio.dataset)

    clinicaldata = getClinicalData(mycgds, cbio.dataset)
    rownames(clinicaldata) = gsub('\\.', '-', rownames(clinicaldata))

    ofile <- paste(cbio.study, cbio.dataset, samples.name, "Rdata", sep = ".")

    ret = NULL
    ret$profile = data
    ret$clinical = clinicaldata

    save(ret, file=ofile)

    cat(paste("\nData exported to file: ", ofile, sep = ""))
    write.table(data, file = ofile)

    return(ret)
}


#### end of file -- loading.R
