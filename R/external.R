#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

#' Create an input file for MUTEX
#' (ref: https://code.google.com/p/mutex/ )
#' @title export,mutex
#'
#' @examples
#' data(gistic)
#' dataset = import.GISTIC(gistic)
#' export.mutex(dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @param filename The name of the file
#' @param filepath The path where to save the file
#' @param label.mutation The event type to use as mutation
#' @param label.amplification The event type to use as amplification (can be a list)
#' @param label.deletion The event type to use as amplification (can be a list)
#' @return A MUTEX example matrix
#' @export export.mutex
#' @importFrom utils write.table
#' 
export.mutex <- function(x, 
                         filename = 'tronco_to_mutex',
                         filepath = './',
                         label.mutation = 'SNV',
                         label.amplification = list('High-level Gain'),
                         label.deletion = list('Homozygous Loss')) {

    is.compliant(x)
    data = x
    alteration =
        list(unlist(label.mutation),
             unlist(label.amplification),
             unlist(label.deletion))

    ## Merge amplification.
    
    if (length(label.amplification) >= 0) {
        amplification = label.amplification[[1]]
    }
    if (length(label.amplification) >= 2) {
        amplification = 'amplification'
        data = join.types(data, label.amplification[[1]], label.amplification[[2]], 'amplification', 'red')
    }
    if (length(label.amplification) > 2) {
        for (label in label.amplification[3:length(label.amplification)]) {
            data = join.types(data, label, 'amplification', 'amplification', 'red')
        }
    }

    ## Merge deletion.
    
    if (length(label.deletion) >= 0) {
        deletion = label.deletion[[1]]
    }
    if (length(label.deletion) >= 2) {
        deletion = 'deletion'
        data = join.types(data, label.deletion[[1]], label.deletion[[2]], 'deletion', 'blue')
    }
    if (length(label.deletion) > 2) {
        for (label in label.deletion[3:length(label.deletion)]) {
            data = join.types(data, label, 'deletion', 'deletion', 'blue')
        }
    }

    ## Merge mutation.
    
    if (length(label.mutation) >= 0) {
        mutation = label.mutation[[1]]
    }
    if (length(label.mutation) >= 2) {
        mutation = 'mutation'
        data = join.types(data, label.mutation[[1]], label.mutation[[2]], 'mutation', 'green')
    }
    if (length(label.mutation) > 2) {
        for (label in label.mutation[3:length(label.mutation)]) {
            data = join.types(data, label, 'mutation', 'mutation', 'green')
        }
    }

    samples = rownames(data$genotypes) 
    genes = unique(data$annotation[,'event'])

    mutex.matrix = matrix(0, nrow = length(genes), ncol = length(samples))
    colnames(mutex.matrix) = samples
    rownames(mutex.matrix) = genes

    ## Legend:
    ## 0: no alteration
    ## 1: mutation
    ## 2: amplification
    ## 4: deletion
    ## 3: 1+2 a+m
    ## 5: 1+4 d+m

    legend = list(1, 2, 4)
    names(legend) = list(mutation, amplification, deletion)
    tronco.matrix = data$genotypes

    for (sample in rownames(tronco.matrix)) {
        for (gene in colnames(tronco.matrix)) {
            type = data$annotations[[gene, 'type']]

            if (type %in% alteration && tronco.matrix[sample, gene] == 1) {
                to.add = legend[[data$annotations[[gene, 'type']]]]
                actual.value = mutex.matrix[data$annotations[[gene, 'event']], sample]
                mutex.matrix[data$annotations[[gene, 'event']], sample] = actual.value + to.add
            }
        }
    }

    ## Reassign value according to mutex notation

    ## Legend:
    ## 0: no alteration
    ## 1: mutation
    ## 2: amplification
    ## 3: deletion
    ## 4: 1+2 a+m
    ## 5: 1+4 d+m

    ## Move a+m to 10
    
    mutex.matrix[which(mutex.matrix == 3)] = 10

    ## move deletion to 3
    
    mutex.matrix[which(mutex.matrix == 4)] = 3

    ## move a+m to 4
    
    mutex.matrix[which(mutex.matrix == 10)] = 4

    mutex.header = append("Symbol", samples)

    filepath = if (grepl("\\/$", filepath)) filepath else paste0(filepath, "/")
    con = paste0(filepath, filename)
    write(mutex.header, file = con, sep = "\t", ncolumns = length(mutex.header))
    write.table(mutex.matrix, con, sep="\t", append = TRUE, col.names = FALSE, quote = FALSE)

    return(mutex.matrix)
}


#' Create a .mat file which can be used with NBS clustering
#' (ref: http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26)
#' @title export.nbs.input
#'
#' @param x A TRONCO compliant dataset.
#' @param map_hugo_entrez Hugo_Symbol-Entrez_Gene_Id map
#' @param file output file name
#' @importFrom R.matlab writeMat
#' @export export.nbs.input
#' 
export.nbs.input <-function(x, 
                            map_hugo_entrez,
                            file = 'tronco_to_nbs.mat') {

    is.compliant(x);

    cat('*** Exporting for NBS v. 0.2\n')
    cat('Preparing binary input matrix\n')

    ## gene_indiv_mat <- the matrix
    gene_indiv_mat = as.matrix(x$genotypes)

    ## Remove colnames and rownames from gene_indiv_mat.
    
    rownames(gene_indiv_mat) = NULL
    colnames(gene_indiv_mat) = NULL

    cat('Preparing samples IDs \n')

    ## sample_id <- patient id
    sample_id = as.samples(x)

    cat('Preparing genes list (should be Hugo_Symbol) \n')
    
    ## gene_id_symbol <- sorted name of events
    gene_id_symbol = as.genes(x)

    cat('Preparing genes map (should be Hugo_Symbol -> Entrez_Gene_Id) \n')
    
    if (!('Hugo_Symbol' %in% colnames(map_hugo_entrez))) {
        stop('No Hugo_Symbol column in the input map: ', colnames(map_hugo_entrez))
    }
    if (!('Entrez_Gene_Id' %in% colnames(map_hugo_entrez))) {
        stop('No Entrez_Gene_Id column in the input map: ', colnames(map_hugo_entrez))
    }

    gene_id_all =
        mapply(function(x) as.numeric(map_hugo_entrez[[which(map_hugo_entrez[ ,'Hugo_Symbol'] == x),
                                                       'Entrez_Gene_Id']]),
               gene_id_symbol)

    file = if (grepl("\\.mat$", file)) file else paste0(file, ".mat")
    con = paste0(file)

    cat('Writing Matlab file to disk:', file,  ' ..... ' )
    writeMat(con,
             gene_indiv_mat = gene_indiv_mat,
             gene_id_all = gene_id_all,
             sample_id = sample_id,
             gene_id_symbol = gene_id_symbol)
    cat('DONE')
}


#' Create a list of unique Mutex groups for a given fdr cutoff
#' current Mutex version is Jan 8, 2015
#' (ref: https://code.google.com/p/mutex/ )
#'
#' @title import.mutex.groups
#' @param file Mutex results ("ranked-groups.txt" file)
#' @param fdr cutoff for fdr
#' @param display print summary table of extracted groups
#' @export import.mutex.groups
#' @importFrom utils count.fields
#' 
import.mutex.groups <- function(file, fdr=.2, display = TRUE) {
    ## Found somewhere on the web - makes sense.
    
    read.irregular <- function(filenm) {
        fileID <- file(filenm,open="rt")
        nFields <- count.fields(fileID)
        mat <- matrix(nrow=length(nFields),ncol=max(nFields))
        invisible(seek(fileID,where=0,origin="start",rw="read"))
        for (i in 1:nrow(mat) ) {
            mat[i, 1:nFields[i]] = scan(fileID, what = "", nlines = 1, quiet = TRUE)
        }
        close(fileID)
        df = data.frame(mat, stringsAsFactors = FALSE)
        return(df)
    }

    x = read.irregular(file)

    ## Check header.
    
    if (any(x[1,1:3] != c('Score', 'q-val', 'Members')))
        warning('File header does not seem to contain \'Score\', \'q-val\' and \'Members field\'.\n',
                'Are you sure this is a Mutex result file?')

    ## Remove header.
    
    cat(paste('*** Groups extracted - ', (nrow(x) -1), ' total groups.\n', sep=''))
    x = x[-1, , drop = FALSE]     # this is c('Score', 'q-val', 'Members')
    x[, 1] = as.numeric(x[,1])          # fdr
    x[, 2] = as.numeric(x[,2])          # q-value

    ## Remove groups  with low fdr.
    
    res = x[which(x[,1] < fdr), , drop = FALSE] 

    ## Remove duplicated groups (permutations).
    
    res.g = res[, 3:ncol(res)]

    for (i in 1:nrow(res.g))
        res[i,3:ncol(res)] = sort(res.g[i,], na.last = TRUE)   
    res = res[!duplicated((res[ , 3:ncol(res), drop = FALSE])), ] 

    cat(paste('Selected ',
              nrow(res),
              ' unique groups with fdr < ',
              fdr,
              '\n',
              sep = ''))

    ## Create groups.
    
    groups <- function(g) {
        g = g[3:length(g)]
        g = g[!is.na(g)]
        names(g) = NULL
        return(sort(g))
    }

    G = list()
    for (i in 1:nrow(res)) {
        gr = list(groups(res[i, ]))
        names(gr) = paste('MUTEX_GROUP', i, sep='')
        G = append(G,gr)  
    }

    rownames(res) = names(G)
    colnames(res)[1:2] = c('fdr', 'score')

    ## Summary report.
    
    if (display) 
        { 
            print(res)
        }
    return(G)
}


#' Check if there are multiple sample in x, according to TCGA barcodes naming
#' @title TCGA.multiple.samples
#' 
#' @examples
#' data(test_dataset)
#' TCGA.multiple.samples(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A list of barcodes. NA if no duplicated barcode is found
#' @export TCGA.multiple.samples
#' 
TCGA.multiple.samples <- function(x) {
    is.compliant(x)

    samples = as.samples(x)
    samples.truncated = substring(samples, 0, 12)

    patients = unique(samples.truncated)

    if (length(patients) != nsamples(x)) {
        dup.samples.start = which(duplicated(samples.truncated)) 
        dup.samples.last = which(duplicated(samples.truncated, fromLast = TRUE))

        return(sort(samples[c(dup.samples.start, dup.samples.last)])) 
    } else
        return(NA)
}


#' If there are multiple sample in x, according to TCGA barcodes naming, remove them
#' @title TCGA.remove.multiple.samples
#' 
#' @examples
#' data(test_dataset)
#' TCGA.remove.multiple.samples(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A TRONCO compliant dataset
#' @export TCGA.remove.multiple.samples
#' 
TCGA.remove.multiple.samples <- function(x) {
    is.compliant(x, err.fun = 'Removing TCGA multiple samples (input)')

    dup = TCGA.multiple.samples(x)
    dup.truncated = substring(dup, 0, 12)
    patients = unique(dup.truncated)

    for (i in 1:length(patients)) {
        patients.samples = which(dup.truncated == patients[i])
        multiple.samples = dup[patients.samples]

        cat('Patient', patients[i], 'with sample aliquotes\n' )
        print(substring(multiple.samples, 14, 29))

        keep = max(multiple.samples)
        discard = multiple.samples[which(multiple.samples != keep)]

        cat('Selecting', keep, '\n')
        x = delete.samples(x, discard)
    }

    is.compliant(x, err.fun = 'Removing TCGA multiple samples (output)')
    return(x)
}


#' Keep only the first 12 character of samples barcode if there are no duplicates
#' @title TCGA.shorten.barcodes
#' 
#' @examples
#' data(test_dataset)
#' TCGA.shorten.barcodes(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A TRONCO compliant dataset
#' @export TCGA.shorten.barcodes
#' 
TCGA.shorten.barcodes <- function(x) {
    is.compliant(x, err.fun='Shartening TCGA barcodes (input)')

    ## Check if it has duplicated barcodes.
    
    if (!all(is.na(TCGA.multiple.samples(x))))
        stop(paste('This dataset contains multiple samples for some patients - cannot consolidate.',
                   '\n Samples with barcodes indicating multiple patients: \n',
                   paste(TCGA.multiple.samples(x), collapse = '\n'),
                   '.',
                   sep = ''))

    ## Shorten sample barcodes.
    
    rownames(x$genotypes) = substring(rownames(x$genotypes), 0, 12)
    if (has.stages(x)) rownames(x$stages) = rownames(x$genotypes)

    is.compliant(x, err.fun='Shartening TCGA barcodes (output)')
    return(x)    
}


#' Map clinical data from the TCGA format
#' @title  TCGA.map.clinical.data
#'
#' @param file A file with the clinical data
#' @param sep file delimiter
#' @param column.samples Required columns
#' @param column.map Map to the required columns
#' @return a map
#' @export TCGA.map.clinical.data
#' @importFrom utils read.delim
#' 
TCGA.map.clinical.data <- function(file, sep='\t', column.samples, column.map) {

    data =
        read.delim(file = file,
                   sep = sep,
                   header = TRUE,
                   stringsAsFactors = FALSE)

    if (!(column.samples %in% colnames(data))) 
        stop(paste('Cannot find samples column \"',
                   column.samples,
                   '\". Available columns: \n\t',
                   paste(colnames(data), collapse = '\n\t'),
                   sep = ''))

    if (!(column.map %in% colnames(data))) 
        stop(paste('Cannot find required map column \"',
                   column.map,
                   '\". Available columns: \n\t',
                   paste(colnames(data), collapse = '\n\t'),
                   sep = ''))

    map = data.frame(data[ , column.map], row.names = data[ , column.samples])
    colnames(map) = column.map

    return(map)
}


## Internal function

sample.RColorBrewer.colors <- function(palette, ncolors) {
    if (!palette %in% rownames(brewer.pal.info))
        stop('Invalid RColorBrewer palette.')

    pmax.cols = brewer.pal.info[palette, 'maxcolors']

    cols = min(pmax.cols , ncolors)
    cols = ifelse(cols < 3, 3, cols)

    colors = brewer.pal(n=cols, name=palette)
    if (ncolors < 3) colors = colors[1:ncolors]
    else colors =  colorRampPalette(colors)(ncolors)

    return(colors)
}

#' Create a graphML object which can be imported in cytoscape
#' This function is based on the tronco.plot fuction
#' 
#' @title  export.graphml
#' 
#' @examples
#' data(test_model)
#' export.graphml(test_model, file='text.xml', scale.nodes=0.3)
#'
#' @param x A TRONCO compliant dataset
#' @param file Where to save the output
#' @param ... parameters for tronco.plot
#' @export export.graphml
#' @importFrom igraph write.graph V V<- set.vertex.attribute
#' @importFrom igraph set.edge.attribute set.graph.attribute
#' @importFrom grDevices rgb col2rgb
#' @importFrom utils packageVersion
#' 
export.graphml <- function(x, file, ...) {

    is.compliant(x)
    is.model(x)

    plot.output = tronco.plot(x, export.igraph = TRUE, ...)
    graph = plot.output$graph
    nodes = plot.output$nodes
    edges = plot.output$edges
    description = plot.output$description
    models = plot.output$models
    node.names = V(graph)$name
    edge.names = names(edges$label)

    ## Display information about vertex name

    V(graph)$name = V(graph)$label

    ## Prepare and save vertex label

    vertex.label = sapply(node.names, function(node) {
        return(gsub("\\", "", nodes$label[node], fixed = TRUE))
    })
    V(graph)$label = vertex.label

    ## Prepare and save vertex type

    vertex.type = sapply(node.names, function(node) {
        if (is.logic.node(node)) {
            return('')
        }
        return(as.events(x)[node,'type'])
    })
    graph = set.vertex.attribute(graph, 'type', value=vertex.type)
    
    ## Prepare and save vertex bg color

    vertex.fillcolor = sapply(node.names, function(node){
        rgb(t(col2rgb(nodes$fillcolor[node])), maxColorValue = 255)
    })
    graph = set.vertex.attribute(graph, 'fillcolor', value=vertex.fillcolor)

    ## Prepare and save vertex font color

    vertex.fontcolor = sapply(node.names, function(node){
        rgb(t(col2rgb(nodes$fontcolor[node])), maxColorValue = 255)
    })
    graph = set.vertex.attribute(graph, 'fontcolor', value=vertex.fontcolor)

    ## Prepare and save vertex border color

    vertex.bordercolor = sapply(node.names, function(node){
        rgb(t(col2rgb(nodes$color[node])), maxColorValue = 255)
    })
    graph = set.vertex.attribute(graph, 'bordercolor', value=vertex.bordercolor)

    ## Prepare and save vertex shape

    vertex.shape = sapply(node.names, function(node){
        shape = nodes$shape[node]
        if (shape == 'box') {
            return('Rectangle')
        }
        return(shape)
    })
    graph = set.vertex.attribute(graph, 'shape', value=vertex.shape)

    ## Prepare and save vertex width

    vertex.width = sapply(node.names, function(node){
        return(nodes$width[node] * 50)
    })
    graph = set.vertex.attribute(graph, 'width', value=vertex.width)

    ## Prepare and save vertex height

    vertex.height = sapply(node.names, function(node){
        return(nodes$height[node] * 50)
    })
    graph = set.vertex.attribute(graph, 'height', value=vertex.height)

    ## Prepare and save label fontsize

    vertex.fontsize = sapply(node.names, function(node){
        return(nodes$fontsize[node] * 1.5)
    })
    graph = set.vertex.attribute(graph, 'fontsize', value=vertex.fontsize)

    ## Prepare and save border width

    vertex.borderwidth = sapply(node.names, function(node){
        return(nodes$lwd[node])
    })
    graph = set.vertex.attribute(graph, 'borderwidth', value=vertex.borderwidth)

    ## Prepare and save vertex shape

    edge.line = sapply(edge.names, function(edge){
        line = edges$lty[edge]
        if (line == 'dashed') {
            return('Dash')
        }
        return('Solid')
    })
    graph = set.edge.attribute(graph, 'line', value=edge.line)

    ## Prepare and save arrow type

    edge.arrow = sapply(edge.names, function(edge){
        arrow = edges$arrowsize[edge]
        if (arrow == 1) {
            return('True')
        }
        return('False')
    })
    graph = set.edge.attribute(graph, 'arrow', value=edge.arrow)

    ## Prepare and save edge color

    edge.color = sapply(edge.names, function(edge){
        rgb(t(col2rgb(edges$color[edge])), maxColorValue = 255)
    })
    graph = set.edge.attribute(graph, 'color', value=edge.color)

    ## Prepare and save graph attributes

    graph = set.graph.attribute(graph, 'name', description)
    graph = set.graph.attribute(graph, 
                                'models',
                                paste(x$parameters$algorithm, models, collapse=' - '))
    graph = set.graph.attribute(graph,
                                'informations',
                                paste0('Generated with TRONCO v', packageVersion('TRONCO')))


    write.graph(graph, file=file, format='graphml')
}


#### end of file -- external.R
