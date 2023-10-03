#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2017, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' oncoPrint : plot a genotype. For details and examples 
#' regarding the visualization through oncoprints, we refer to the Vignette Section 4.4. 
#' 
#' @title oncoprint
#' @param x A TRONCO compliant dataset
#' @param excl.sort Boolean value, if TRUE sorts samples to enhance exclusivity of alterations
#' @param samples.cluster Boolean value, if TRUE clusters samples (columns). Default FALSE
#' @param genes.cluster Boolean value, if TRUE clusters genes (rows). Default FALSE
#' @param file If not NA write to \code{file} the Oncoprint, default is NA (just visualization).
#' @param ann.stage Boolean value to annotate stage classification, default depends on \code{x}
#' @param ann.hits Boolean value to annotate the number of events in each sample, default is TRUE
#' @param stage.color RColorbrewer palette to color stage annotations. Default is 'YlOrRd'
#' @param hits.color RColorbrewer palette to color hits annotations. Default is 'Purples'
#' @param null.color Color for the Oncoprint cells with 0s, default is 'lightgray'
#' @param border.color Border color for the Oncoprint, default is white' (no border)
#' @param text.cex Title and annotations cex, multiplied by font size 7
#' @param font.column If NA, half of font.row is used
#' @param font.row If NA, max(c(15 * exp(-0.02 * nrow(data)), 2)) is used, where data is the data 
#' visualized in the Oncoprint
#' @param title Oncoprint title, default is as.name(x) - see \code{as.name}
#' @param sample.id If TRUE shows samples name (columns). Default is FALSE
#' @param hide.zeroes If TRUE trims data - see \code{trim} - before plot. Default is FALSE 
#' @param legend If TRUE shows a legend for the types of events visualized. Defualt is TRUE
#' @param legend.cex Default 0.5; determines legend size if \code{legend = TRUE}
#' @param cellwidth Default NA, sets autoscale cell width 
#' @param cellheight Default NA, sets autoscale cell height
#' @param group.by.label Sort samples (rows) by event label - usefull when multiple events per gene are
#' available 
#' @param group.samples If this samples -> group map is provided, samples are grouped as of groups
#' and sorted according to the number of mutations per sample - usefull when \code{data} was clustered
#' @param group.by.stage Default FALSE; sort samples by stage.
#' @param gene.annot Genes'groups, e.g. list(RAF=c('KRAS','NRAS'), Wnt=c('APC', 'CTNNB1')). Default is NA.
#' @param gene.annot.color Either a RColorColorbrewer palette name or a set of custom colors matching names(gene.annot)
#' @param show.patterns If TRUE shows also a separate oncoprint for each pattern. Default is FALSE
#' @param annotate.consolidate.events Default is FALSE. If TRUE an annotation for events to consolidate is shown.
#' @param txt.stats By default, shows a summary statistics for shown data (n,m, |G| and |P|)
#' @param gtable If TRUE return the gtable object
#' @param ... other arguments to pass to pheatmap
#' @export oncoprint
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom gtable gtable gtable_add_grob gtable_height gtable_width
#' @importFrom grDevices colorRampPalette
#' 
oncoprint <- function(x, 
                      excl.sort = TRUE, 
                      samples.cluster = FALSE, 
                      genes.cluster = FALSE, 
                      file = NA, 
                      ann.stage = has.stages(x), 
                      ann.hits = TRUE, 
                      stage.color = 'YlOrRd', 
                      hits.color = 'Purples',  
                      null.color = 'lightgray', 
                      border.color = 'white', 
                      text.cex = 1.0, 
                      font.column = NA, 
                      font.row = NA, 
                      title = as.description(x),
                      sample.id = FALSE,
                      hide.zeroes = FALSE,
                      legend = TRUE,
                      legend.cex = 0.5,
                      cellwidth = NA, 
                      cellheight = NA,
                      group.by.label = FALSE,
                      group.by.stage = FALSE,
                      group.samples = NA,
                      gene.annot = NA,
                      gene.annot.color = 'Set1',
                      show.patterns = FALSE,
                      annotate.consolidate.events = FALSE,
                      txt.stats =
                          paste(nsamples(x),
                                ' samples\n',
                                nevents(x),
                                ' events\n',
                                ngenes(x),
                                ' genes\n',
                                npatterns(x),
                                ' patterns',
                                sep = ''),
                      gtable = FALSE,
                      ...) {

    font.size = text.cex * 7

    ##  This function sorts a matrix to enhance mutual exclusivity
    
    exclusivity.sort <- function(M) {
        geneOrder <- sort(rowSums(M),
                          decreasing = TRUE,
                          index.return = TRUE)$ix
        scoreCol <- function(x) {
            score <- 0;
            for (i in 1:length(x)) {
                if(x[i]) {
                    score <- score + 2^(length(x)-i);
                }
            }
            return(score);
        }
        scores <- apply(M[geneOrder, , drop = FALSE ], 2, scoreCol);
        sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix
        res = list()
        res$geneOrder = geneOrder
        res$sampleOrder = sampleOrder
        res$M = M[geneOrder, sampleOrder]

        return(res);
    }

    ## Check input data.
    
    cat(paste('*** Oncoprint for "',
              title,
              '"\nwith attributes: stage = ',
              ann.stage,
              ', hits = ',
              ann.hits,
              '\n',
              sep = ''))
    is.compliant(x, 'oncoprint', stage=ann.stage)
    x = enforce.numeric(x)

    ## If hide.zeros trim x
    
    if (hide.zeroes) {
        cat(paste('Trimming the input dataset (hide.zeroes).\n',
                  sep = ''))
        x = trim(x)
    }

    ##  Reverse the heatmap under the assumption that ncol(data) << nrow(data)
    data = t(x$genotypes)
    nc = ncol(data) 
    nr = nrow(data)

    ## Sort data, if required. excl.sort and group.samples are not compatible
    hasGroups = !any(is.na(group.samples))

    if (group.by.stage && !ann.stage)
        stop('Cannot group samples by stage if no stage annotation is provided.')

    if (group.by.stage && excl.sort)
        if (excl.sort && (hasGroups || group.by.stage)) 
            stop('Disable sorting for mutual exclusivity (excl.sort=FALSE) or avoid using grouped samples.')

    ## Exclusivity sort
    if (excl.sort && nevents(x) > 1) {
        cat(paste('Sorting samples ordering to enhance exclusivity patterns.\n',
                  sep = ''))
        sorted.data = exclusivity.sort(data)
        data = sorted.data$M  
    }

    if (group.by.stage) {  
        ord.stages = as.stages(x)[order(as.stages(x)[,1]), , drop = FALSE]
        cat('Grouping samples by stage annotation.\n')

        aux.fun <- function(samp) {
            print(samp)
            sub.data = data[, samp, drop= FALSE]      
            sub.data = sub.data[, order(colSums(sub.data), decreasing = FALSE), drop = FALSE]
            return(sub.data)
        }

        new.data = NULL
        u.stages = sort(unlist(unique(as.stages(x))), na.last = TRUE)

        ## print(str(ord.stages))

        for (i in u.stages) {
            print(ord.stages[which(ord.stages == i), , drop=FALSE])
            new.data = cbind(new.data, aux.fun(rownames(ord.stages[which(ord.stages == i), , drop= FALSE])))
        }
        data = new.data
        data = data[order(rowSums(data), decreasing = TRUE), , drop = FALSE ]
    }

    ## Samples grouping via hasGroups
    if (hasGroups) {
        group.samples[, 1] = as.character(group.samples[,1])
        grn = rownames(group.samples)

        cat(paste('Grouping samples according to input groups (group.samples).\n', sep = ''))
        if (any(is.null(grn)))
            stop('"group.samples" should be matrix with sample names and group assignment.')

        if (!setequal(grn, as.samples(x)))
            stop(paste0('Missing group assignment for samples: ',
                        paste(setdiff(as.samples(x), grn),
                              collapse=', '),
                        '.'))

        ## Order groups by label, and then data (by column)
        order = order(group.samples[,1])
        group.samples = group.samples[order, , drop = FALSE]

        data = data[, rownames(group.samples)]            
        data = data[order(rowSums(data), decreasing = TRUE), , drop = FALSE]    

        groups = unique(group.samples[,1])

        for(i in 1:length(groups)) {
            subdata = data[, group.samples == groups[i], drop = FALSE]
            subdata = subdata[, order(colSums(subdata), decreasing = TRUE), drop = FALSE]
            data[ , group.samples == groups[i]] = subdata
        }
    } 

    ## If group.by.label group events involving the gene symbol
    if (group.by.label) {
        cat(paste('Grouping events by gene label, samples will not be sorted.\n',
                  sep = ''))
        genes = as.genes(x)
        data = data[ order(x$annotations[rownames(data), 'event']), ]
    }

    cn = colnames(data)
    rn = rownames(data)

    ## SAMPLES annotations: hits (total 1s per sample), stage or groups
    samples.annotation = data.frame(row.names = cn,  stringsAsFactors = FALSE)
    nmut = colSums(data)

    if (ann.hits)  samples.annotation$hits = nmut
    if (ann.stage) samples.annotation$stage = as.stages(x)[cn, 1]
    if (hasGroups) samples.annotation$cluster = group.samples[cn, 1]

    ## Color samples annotation
    annotation_colors = NULL

    ## Color hits
    if (ann.hits) {
        hits.gradient = (colorRampPalette(brewer.pal(6, hits.color))) (max(nmut))
        annotation_colors = append(annotation_colors, list(hits=hits.gradient))
    }

    ## Color stage
    if (ann.stage) { 
        cat('Annotating stages with RColorBrewer color palette',
            stage.color,
            '\n')

        different.stages = sort(unique(samples.annotation$stage))

        stage.color.attr = sample.RColorBrewer.colors(stage.color, length(different.stages))    
        stage.color.attr = append(stage.color.attr, "#FFFFFF")      

        samples.annotation[is.na(samples.annotation)] = "none"
        names(stage.color.attr) = append(different.stages, "none")    
        annotation_colors = append(annotation_colors, list(stage=stage.color.attr))
    }

    ## Color groups
    if (hasGroups) {
        ngroups = length(unique(group.samples[,1]))
        cat('Grouping labels:', paste(unique(group.samples[,1]), collapse=', '), '\n')

        group.color.attr = sample.RColorBrewer.colors('Accent', ngroups)
        names(group.color.attr) = unique(group.samples[,1])
        annotation_colors = append(annotation_colors, list(cluster=group.color.attr))
    }

    ## GENES/EVENTS annotations: groups or indistinguishable.
    genes.annotation = NA

    ## Annotate genes groups.
    
    if (!all(is.na(gene.annot))) {
        names = names(gene.annot)

        genes.annotation = data.frame(row.names = rn, stringsAsFactors = FALSE)
        genes.annotation$group = rep("none", nrow(data))

        for (i in 1:length(names)) {
            pathway = names[i]
            genes.pathway = rownames(as.events(x, genes=gene.annot[[names[i]]]))
            genes.annotation[genes.pathway, 'group'] = names[i] 
        }

        if (length(gene.annot.color) == 1
            && gene.annot.color %in% rownames(brewer.pal.info)) {
            cat('Annotating genes with RColorBrewer color palette', gene.annot.color, '.\n')
            gene.annot.color = sample.RColorBrewer.colors(gene.annot.color, length(names))
            gene.annot.color = append(gene.annot.color, "#FFFFFF")      
        } else {
            if (length(gene.annot.color) != length(names)) 
                stop('You did not provide enough colors to annotate',
                     length(names),
                     'genes.\n',
                     'Either set gene.annot.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')

            cat('Annotating pathways with custom colors:',
                paste(gene.annot.color, collapse=', '),
                '\n')
            gene.annot.color = append(gene.annot.color, "#FFFFFF")
        }
        names(gene.annot.color) = append(names, "none")
        gene.annot.color = gene.annot.color[ unique(genes.annotation$group) ]   
        annotation_colors = append(annotation_colors, list(group=gene.annot.color))
    }   

    ## Annotate events to consolidate
    if (annotate.consolidate.events) {
        cat('Annotating events to consolidate - see consolidate.data\n')
        invalid = consolidate.data(x, print = FALSE)

        genes.annotation$consolidate = 'none'

        if (length(invalid$indistinguishable) > 0) {

            genes.annotation[
                             unlist(
                                 lapply(
                                     invalid$indistinguishable, 
                                     function(z) {
                                         return(rownames(unlist(z)))
                                     }
                                 )),
                             'consolidate'] = 'indistinguishable'
        }

        if (length(invalid$zeroes) > 0) {
            genes.annotation[ unlist(invalid$zeroes), 'consolidate'] = 'missing'
        }

        if (length(invalid$ones) > 0) {
            genes.annotation[ unlist(invalid$ones), 'consolidate'] = 'persistent'
        }

        consolidate.colors = c('white', 'brown1', 'darkorange4', 'firebrick4', 'deepskyblue3')
        names(consolidate.colors) = c('none', 'indistinguishable', 'missing', 'persistent')
        consolidate.colors = consolidate.colors[unique(genes.annotation$consolidate)]

        annotation_colors =
            append(annotation_colors,
                   list(consolidate = consolidate.colors)
                   )
    }

    ## Augment gene names with frequencies
    genes.freq = rowSums(data) / nsamples(x)
    gene.names = x$annotations[rownames(data), 2]
    gene.names =
        paste(round(100 * genes.freq, 0),
              '% ',
              gene.names,
              sep = '') # row labels


    ## Augment data to make type-dependent colored plots.
    data.lifting <- function(obj, matrix) {
        types = as.types(obj)
        map.gradient = null.color

        for (i in 1:ntypes(obj)) {
            events = as.events(obj, types = as.types(obj)[i])
            keys = rownames(events)

            if (ntypes(obj) > 1) {
                keys.subset =
                    keys[unlist(lapply(keys,
                                       function(obj, matrix) {

                                           ## Are you sure (obj %in% # matrix)
                                           ## is not sufficient?

                                           if (obj %in% matrix) {
                                               TRUE
                                           } else {
                                               FALSE
                                           }},
                                       rownames(matrix)))]
                sub.data = matrix[keys.subset, , drop = FALSE]

                ## shift 1s to 'i', add color to the map.
                
                idx = which(sub.data == 1)
                if(length(idx) > 0) map.gradient = cbind(map.gradient, as.colors(obj)[i])

                sub.data[idx] = i
                matrix[keys.subset, ] = sub.data 

            } else {
                map.gradient = cbind(map.gradient, as.colors(obj)[i])
            }
        }

        map.gradient = c(null.color, as.colors(x))
        names(map.gradient)[1] = 'none'

        return(list(data=matrix, colors=map.gradient))
    }

    pheat.matrix = data.lifting(x,data)
    map.gradient = pheat.matrix$colors
    data = pheat.matrix$data

    ## Set fontisize col/row.
    if (is.na(font.row)) {
        font.row = max(c(15 * exp(-0.02 * nrow(data)), 2))    
        cat(paste('Setting automatic row font (exponential scaling): ',
                  round(font.row, 1),
                  '\n',
                  sep = ''))
    }
    
    if (is.na(font.column) && sample.id) {
        font.column = font.row/2    
        cat(paste('Setting automatic samples font half of row font: ',
                  round(font.column, 1),
                  '\n',
                  sep = '')) 
    }

    ## Finalizing legends etc.
    legend.labels = c('none', unique(x$annotations[, 1]))    
    legend.labels = legend.labels[1:(max(data) + 1)]

    if (samples.cluster) cat('Clustering samples and showing dendogram.\n')
    if (genes.cluster) cat('Clustering alterations and showing dendogram.\n')

    if (is.null(annotation_colors)) annotation_colors = NA

    if (length(list(...)) > 0) {
        cat('Passing the following parameters to TRONCO pheatmap:\n')
        print(list(...))
    }

    ## Real pheatmap.
    
    if (ncol(samples.annotation) == 0) {
        samples.annotation = NA
    }

    ret =
        pheatmap(data, 
                 scale = "none", 
                 color = map.gradient, 
                 cluster_cols = samples.cluster,
                 cluster_rows = genes.cluster,
                 main = title,
                 fontsize = font.size,
                 fontsize_col = font.column,
                 fontsize_row = font.row,
                 annotation_col = samples.annotation,
                 annotation_row = genes.annotation,
                 annotation_colors = annotation_colors, 
                 border_color = border.color,
                 border = TRUE,
                 margins = c(10,10),
                 cellwidth = cellwidth, 
                 cellheight = cellheight,
                 legend = legend,             
                 legend_breaks = c(0:max(data)),
                 legend_labels = legend.labels,
                 legend.cex = legend.cex,
                 labels_row = gene.names,
                 drop_levels = TRUE,
                 show_colnames = sample.id,
                 filename = file,
                 txt.stats = txt.stats,
                 ...
                 )

    ## Extra patterns   

    patt.table =
        gtable(widths = unit(c(7, 2), "null"),
               heights = unit(c(2, 7), "null"))

    patt.table = list()
    map.gradient = map.gradient[names(map.gradient) != 'Pattern']

    ## print(data)
    
    if (npatterns(x) > 0 && show.patterns) {
        cat('Plotting also', npatterns(x), 'patterns\n')
        patterns = as.patterns(x)

        for (i in 1:length(patterns)) {
            genes.patt = as.events.in.patterns(x, patterns[i]) 
            genes.patt.genos = data[rownames(genes.patt), , drop = FALSE]


            genes.patt.genos.gtable =
                pheatmap(exclusivity.sort(genes.patt.genos)$M,
                         scale = "none",
                         color = map.gradient,
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         main = paste('Pattern:', patterns[i]),
                         fontsize = font.size * .75,
                         fontsize_col = font.column * .75,
                         fontsize_row = font.row * .75,
                         ## annotation_row = genes.annotation,
                         annotation_colors = annotation_colors,
                         annotation_legend  = FALSE,
                         border_color = border.color,
                         border = TRUE,
                         cellwidth = 6,
                         cellheight = 6,
                         legend = FALSE,
                         labels_row = genes.patt[rownames(genes.patt.genos), 'event'],
                         drop_levels = TRUE,
                         show_colnames = FALSE,
                         silent = TRUE,
                         ...)$gtable

            ## patt.table = gtable_add_grob(patt.table, genes.patt.genos.gtable, 2, 1)
            patt.table = append(patt.table, list(genes.patt.genos.gtable))
        }

        str =
            paste('grid.arrange(ret$gtable, arrangeGrob(', 
                  paste('patt.table[[', 1:length(patt.table), ']],', collapse = ''),
                  'ncol = 1), ncol = 1, heights = c(4,',
                  paste(rep(1, length(patt.table) ), collapse=', '),
                  '))')

        str =
            paste('grid.arrange(ret$gtable, arrangeGrob(', 
                  paste('patt.table[[', 1:length(patt.table), ']],', collapse = ''),
                  'ncol = 1), ncol = 1, heights = c(',
                  paste(rep(1, length(patt.table) + 1 ),
                        collapse=', '),
                  '))')

        eval(parse(text = str))
    }


    ##  print(ncol(ret$gtable))
    ##  print(ncol(patt.table))
    
    ## ret$gtable = rbind(ret$gtable, patt.table)
    ## ret$gtable = gtable_add_grob(ret$gtable, patt.table, 1 , 1)
    ## grid.newpage()
    ## grid.draw(ret$gtable)
    if (gtable) {
        return(ret)
    }
}


##### Pathway print
#' Visualise pathways informations
#' @title pathway.visualization
#'
#' @param x A TRONCO complian dataset
#' @param title Plot title
#' @param file To generate a PDF a filename have to be given
#' @param pathways.color A RColorBrewer color palette
#' @param aggregate.pathways Boolean parameter
#' @param pathways Pathways
#' @param ... Additional parameters
#' @return plot information
#' @export pathway.visualization
#' 
pathway.visualization <- function(x, 
                                  title =
                                      paste('Pathways:',
                                            paste(names(pathways), collapse=', ', sep='')),
                                  file = NA,
                                  pathways.color = 'Set2',
                                  aggregate.pathways,
                                  pathways,
                                  ...) { 
    names = names(pathways)

    if (length(pathways.color) == 1
        && pathways.color %in% rownames(brewer.pal.info)) {
        cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
        n = length(names)
        if (n < 3) {
            n = 3
        }
        pathway.colors = brewer.pal(n = n, name=pathways.color)
    } else { 
        if (length(pathways.color) != length(names)) 
            stop('You did not provide enough colors to annotate ',
                 length(names),
                 ' pathways.\n',
                 'Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')

        cat('Annotating pathways with custom colors', paste(pathways.color, collapse=','), '.\n')
        pathway.colors = pathways.color
    }

    names(pathway.colors) = names

    cat(paste('*** Processing pathways: ',
              paste(names, collapse = ', ', sep = ''),
              '\n',
              sep = ''))

    cat(paste('\n[PATHWAY \"', names[1],'\"] ',
              paste(pathways[[1]], collapse = ', ', sep = ''),
              '\n',
              sep = ''))

    data.pathways =
        as.pathway(x,
                   pathway.genes = pathways[[1]], 
                   pathway.name = names[1],
                   aggregate.pathway = aggregate.pathways)

    data.pathways = change.color(data.pathways, 'Pathway', pathway.colors[1])
    data.pathways = rename.type(data.pathways, 'Pathway', names[1])

    if (length(names) > 1) {  
        for(i in 2:length(pathways)) {
            cat(paste('\n\n[PATHWAY \"', names[i],'\"] ', paste(pathways[[i]], collapse=', ', sep=''), '\n', sep=''))
            pathway = as.pathway(x, pathway.genes=pathways[[i]], pathway.name=names[i], aggregate.pathway = aggregate.pathways)

            pathway = change.color(pathway, 'Pathway', pathway.colors[i])
            pathway = rename.type(pathway, 'Pathway', names[i])
            data.pathways = ebind(data.pathways, pathway)
        }
    }

    ret = oncoprint(trim(data.pathways), title=title, file=file, ...)

    return(ret)
}


#' export input for cbio visualization at http://www.cbioportal.org/public-portal/oncoprinter.jsp
#' @title oncoprint.cbio
#'
#' @examples
#' data(crc_gistic)
#' gistic = import.GISTIC(crc_gistic)
#' oncoprint.cbio(gistic)
#'
#' @param x A TRONCO compliant dataset.
#' @param file name of the file where to save the output
#' @param hom.del type of Homozygous Deletion
#' @param het.loss type of Heterozygous Loss
#' @param gain type of Gain
#' @param amp type of Amplification
#' @return A file containing instruction for the CBio visualization Tool
#' @export oncoprint.cbio
#' 
oncoprint.cbio <- function(x, 
                           file = 'oncoprint-cbio.txt', 
                           hom.del = 'Homozygous Loss',
                           het.loss = 'Heterozygous Loss', 
                           gain = 'Low-level Gain', 
                           amp = 'High-level Gain') {
    is.compliant(x)

    ## r = paste(paste(rownames(x$genotypes)), x$annotations[x$annotations[,''],], 'xxx')

    r = 'Sample\tGene\tAlteration\n'
    for (i in 1:nrow(x$genotypes)) {
        for (j in 1:ncol(x$genotypes)) {
            if (x$genotypes[i,j] == 1)
                {
                    s = rownames(x$genotypes)[i]
                    g = x$annotations[colnames(x$genotypes)[j], 'event']

                    t = x$annotations[colnames(x$genotypes)[j], 'type']

                    t.o = 'xxx'
                    if (t == hom.del) t.o = 'HOMDEL'
                    if (t == het.loss)  t.o = 'HETLOSS'
                    if (t == gain)  t.o = 'GAIN'
                    if (t == amp)  t.o = 'AMP'

                    ## cat(paste( s,  g,  t.o, '\n', sep=' ', collpase=''))
                    r = paste(r, s, '\t', g, '\t', t.o, '\n', sep='', collpase='')
                }
        }
    }
    
    write(r, file)
}


#' Generate PDF and laex tables
#' @title genes.table.report
#'
#' @param x A TRONCO compliant dataset.
#' @param name filename
#' @param dir working directory
#' @param maxrow maximum number of row per page
#' @param font document fontsize
#' @param height table height
#' @param width table width
#' @param fill fill color
#' @param silent A parameter to disable/enable verbose messages.
#' @return LaTEX code
#' @importFrom gridExtra grid.table
#' @importFrom xtable xtable
#' @importFrom grDevices pdf dev.cur dev.off dev.set
#' @export genes.table.report
#' 
genes.table.report <- function(x,
                               name,
                               dir = getwd(),
                               maxrow = 33, 
                               font = 10,
                               height = 11,
                               width = 8.5,
                               fill = "lightblue",
                               silent = FALSE) {
    ## Print table with gridExtra and xtables.
    
    print.table <- function(table,
                            name,
                            dir = getwd(),
                            maxrow,
                            font,
                            height, width,
                            fill) {
        cat('\nPrinting PDF and Latex table to files: \n')
        cat(paste('PDF \t\t', dir, '/', name, '.genes-table.pdf\n', sep=''))
        cat(paste('Latex\t\t', dir, '/', name, '.genes-table.tex\n', sep=''))
        cat('\n')

        ## Output pdf
        ## require(gridExtra)  
        ## require(xtable)  

        cur.dev = dev.cur()

        pdf(file = paste(dir, '/', name, '.genes-table.pdf', sep = ''),
            height = height,
            width = width)

        ## Max rows per page
        
        npages = ceiling(nrow(table)/maxrow); 

        #flush.console()

        #pb = txtProgressBar(1, npages, style = 3);      
        for (i in 1:npages) {
            if (!silent) {
              cat('.')
            }
            #setTxtProgressBar(pb, i)  
            idx = seq(1+((i-1)*maxrow), i*maxrow); 

            if (max(idx) > nrow(table)) idx = idx[idx < nrow(table)]   

            grid.newpage()
            grid.table(table[idx, ],
                       gpar.coretext = gpar(fontsize = font),
                       gpar.corefill = gpar(fill = fill, alpha=0.5, col = NA),
                       h.even.alpha = 0.5)
        } 
        #close(pb)
        if (!silent) {
          cat('\n')
        }

        ## Output latex.
        
        print(xtable(table, digits = 0),
              file = paste(dir, '/', name, '.genes-table.tex', sep = ''),
              type = 'latex')

        dev.off()
        dev.set(which = cur.dev)
    }

    cat(paste('Preparing output table with ', ngenes(x), ' genes ...\n'))
    genes = as.genes(x)
    types = as.types(x)

    data = matrix(0, nrow = ngenes(x), ncol = ntypes(x))

    genes.table = data.frame(data, row.names = genes, stringsAsFactors = FALSE)
    colnames(genes.table) = types

    x = enforce.numeric(x)

    #pb = txtProgressBar(1, ngenes(x), style = 3);
    for (i in 1:ngenes(x)) {
        #setTxtProgressBar(pb, i)  
        if (!silent) {
            cat('.')
        }

        g = as.gene(x, genes=genes[i])

        if (ncol(g) > 0) {
            gg = colSums(apply(g, 2, as.numeric))

            genes.table[rownames(genes.table)[i], colnames(g)] = gg
            genes.table[rownames(genes.table)[i], 'Alterations'] =
                paste( round(sum(gg) / nsamples(x) * 100), '%', sep='')  
            genes.table[rownames(genes.table)[i], 'Frequency'] =
                sum(gg) / nsamples(x)  
        }
    }
                                        ## Close progress bar.
    #close(pb)

    if (!silent) {
        cat('\n')
    }

    genes.table = genes.table[order(genes.table$Frequency, decreasing = TRUE), ]
    genes.table$Frequency = NULL

    print.table(table = genes.table,
                name = name,
                dir = getwd(),
                maxrow = maxrow,
                font = font,
                height = height,
                width = width,
                fill = fill)

    return(genes.table)
}


#' @importFrom RColorBrewer brewer.pal
cluster.sensitivity <- function(cluster.map,
                                  reference,
                                  stages = NA,
                                  file = NA) {

    if (ncol(cluster.map) == 1)
        stop('No clustering stability for a unique clustering map!')


    if (!reference %in% colnames(cluster.map)) 
        stop(paste0('The reference cluster specified is not any of: ', 
                    paste(colnames(cluster.map), collapse = ', '), '.'))

    ref.clust = which(reference == colnames(cluster.map))  
    colnames(cluster.map)[ref.clust] =
        paste(colnames(cluster.map)[ref.clust],' [reference]',sep = '') 

    ## Transpose data.
    ## cluster.map = t(cluster.map)

    ## Sort data according to reference row.
    
    cluster.map =
        cluster.map[sort(cluster.map[, ref.clust ],
                         decreasing = FALSE,
                         index.return = TRUE)$ix, ];

    ## Get unique clustering IDs.
    
    id = apply(cluster.map, 2, unique)
    id = unique(unlist(id))

    print(apply(cluster.map, 2, unique))

    ## Compute the clustering score.
    
    subdata = cluster.map[, -ref.clust]
    refcol = cluster.map[, ref.clust]
    urefcol = unique(refcol)

    cat('Found the following cluster labels:', urefcol, '\n')

    cat('Computing clustering scores ... ')

    score = rep(0, nrow(cluster.map))

    for (i in 1:length(urefcol)) {
        tmp = as.matrix(subdata[which(refcol == i), ]);

        curr.score = 0;
        for (j in 1:ncol(tmp)) {
            curr.cardinality = sort(table(tmp[, j]), decreasing = TRUE)[[1]];
            curr.score =
                curr.score + (nrow(tmp) - curr.cardinality) / nrow(tmp);     
        }
        score[which(refcol == i)] = 1 - (curr.score/nrow(tmp))
    }
    cat('DONE\n')


    ## Create annotations.
    
    cn = rownames(cluster.map)

    annotation =
        data.frame(sensitivity = score,
                   row.names = cn,
                   stringsAsFactors = FALSE)

    if (!all(is.na(stages)))
        annotation$stage = stages[cn,1]

    ## Create colors.
    
    col = brewer.pal(n = length(id), name = 'Set1')

    different.stages = sort(unique(annotation$stage))
    num.stages = length(different.stages)

    stage.color = append(brewer.pal(n = num.stages, name = 'YlOrRd'), '#FFFFFF')
    names(stage.color) = append(levels(as.factor(different.stages)), NA)

    score.color = brewer.pal(n = 3, name = 'Greys')

    ## Annotation colors.
    
    annotation_colors = list(stage=stage.color, sensitivity=score.color)

    ## Settings.
    
    main =
        paste0("Clustering Sensitivity\n Reference : ", 
               reference,
               '\nAgainst : ',
               paste(colnames(subdata), collapse = ', ')) 


    cat('Clustering rows in', nrow(cluster.map), 'clusters.\n')
    order =
        sort(cluster.map[,ref.clust],
             decreasing = FALSE,
             index.return = TRUE) 



    pheatmap(cluster.map,
             scale = "none",
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             color = col,
             main = main,
             fontsize = 6,
             fontsize_col = 8,
             fontsize_row = 2,
             annotation_row = annotation,
             annotation_colors = annotation_colors, 
             border_color = 'lightgray',
             border = TRUE,
             margins = c(10, 10),
             cellwidth = 25,
             cellheight = 2.2,
             ## legend = TRUE,
             legend_breaks = 1:4,
             filename = file,
             gaps_col = 1:3,
             ## display_numbers = T,
             cutree_rows = nrow(cluster.map),
             gaps_row = (match(urefcol, order$x) - 1)
             ## number_format = '%d',                 
             )

    return(cluster.map)
}

lo <- function(rown,
               coln,
               nrow,
               ncol,
               cellheight = NA,
               cellwidth = NA,
               treeheight_col,
               treeheight_row,
               legend,
               annotation_row,
               annotation_col,
               annotation_colors,
               annotation_legend,
               main,
               fontsize,
               fontsize_row,
               fontsize_col,
               gaps_row,
               gaps_col,
               ...) {
    
    ## Get height of colnames and length of rownames.

    if (!is.null(coln[1])) {
        t = c(coln, colnames(annotation_row))
        longest_coln = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_col, ...)
        coln_height =
            unit(1,
                 "grobheight",
                 textGrob(t[longest_coln],
                          rot = 90,
                          gp = do.call(gpar, gp))) + unit(10, "bigpts")
    } else {
        coln_height = unit(5, "bigpts")
    }

    if (!is.null(rown[1])) {
        t = c(rown, colnames(annotation_col))
        longest_rown = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_row, ...)
        rown_width =
            unit(1,
                 "grobwidth",
                 textGrob(t[longest_rown],
                          gp = do.call(gpar, gp))) + unit(10, "bigpts")
    } else {
        rown_width = unit(5, "bigpts")
    }

    gp = list(fontsize = fontsize, ...)
    
    ## Legend position

    if (!is.na(legend[1])) {
        longest_break = which.max(nchar(names(legend)))
        longest_break =
            unit(1.1,
                 "grobwidth",
                 textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
        title_length =
            unit(1.1,
                 "grobwidth",
                 textGrob("Scale", gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    } else {
        legend_width = unit(0, "bigpts")
    }

    ## Set main title height.
    
    if (is.na(main)){
        main_height = unit(0, "npc")
    } else {
        main_height =
            unit(1.5,
                 "grobheight",
                 textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }

    ## Column annotations.
    
    textheight = unit(fontsize, "bigpts")

    if (!is.na(annotation_col[[1]][1])) {
        ## Column annotation height.
        annot_col_height =
            ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

        ## Width of the correponding legend.
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col)) 
        annot_col_legend_width =
            unit(1.2,
                 "grobwidth",
                 textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_col_legend_width = unit(0, "npc")
        }
    } else {
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }

    ## Row annotations.
    if (!is.na(annotation_row[[1]][1])) {
        
        ## Row annotation width.
        annot_row_width =
            ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

        ## Width of the correponding legend.
        
        t =
            c(as.vector(as.matrix(annotation_row)),
              colnames(annotation_row)) 
        annot_row_legend_width =
            unit(1.2,
                 "grobwidth",
                 textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        
        if (!annotation_legend) {
            annot_row_legend_width = unit(0, "npc")
        }
    } else {
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }

    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)

    ## Tree height.
    
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 

    ## Set cell sizes.
    
    if (is.na(cellwidth)) {
        mat_width =
            unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width 
    } else {
        mat_width =
            unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(4, "bigpts")
    }

    if (is.na(cellheight)) {
        mat_height = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
    } else {
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(4, "bigpts")
    }    

    ## Produce gtable.
    
    gt =
        gtable(widths =
                   unit.c(treeheight_row,
                          annot_row_width, 
                          mat_width,
                          rown_width,
                          legend_width,
                          annot_legend_width),
               heights =
                   unit.c(main_height,
                          treeheight_col,
                          annot_col_height,
                          mat_height,
                          coln_height),
               vp = viewport(gp = do.call(gpar, gp)))

    ## print(unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height))

    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, "bigpts")), "bigpts", valueOnly = TRUE) / ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, "bigpts")), "bigpts", valueOnly = TRUE) / nrow

    ## Return minimal cell dimension in bigpts to decide if borders
    ## are drawn.
    
    mindim = min(cw, ch)

    res = list(gt = gt, mindim = mindim)

    return(res)
}


find_coordinates <- function(n, gaps, m = 1:n) {
    if (length(gaps) == 0) {
        return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
    }

    if (max(gaps) > n) {
        stop("Gaps do not match with matrix size")
    }

    size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))

    gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum) 
    coord = m * size + (gaps2 * unit("4", "bigpts"))

    return(list(coord = coord, size = size))
}


draw_dendrogram <- function(hc, gaps, horizontal = TRUE) {
    h = hc$height / max(hc$height) / 1.05
    m = hc$merge
    o = hc$order
    n = length(o)

    m[m > 0] = n + m[m > 0] 
    m[m < 0] = abs(m[m < 0])

    dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
    dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)

    for (i in 1:nrow(m)) {
        dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
        dist[n + i, 2] = h[i]
    }

    draw_connection <- function(x1, x2, y1, y2, y) {
        res =
            list(x = c(x1, x1, x2, x2),
                 y = c(y1, y, y, y2)
                 )
        return(res)
    }

    x = rep(NA, nrow(m) * 4)
    y = rep(NA, nrow(m) * 4)
    id = rep(1:nrow(m), rep(4, nrow(m)))

    for (i in 1:nrow(m)) {
        c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
        k = (i - 1) * 4 + 1
        x[k : (k + 3)] = c$x
        y[k : (k + 3)] = c$y
    }

    x = find_coordinates(n, gaps, x * n)$coord
    y = unit(y, "npc")

    if (!horizontal) {
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
    }
    res = polylineGrob(x = x, y = y, id = id)

    return(res)
}


draw_matrix <- function(matrix,
                        border_color,
                        gaps_rows,
                        gaps_cols,
                        fmat,
                        fontsize_number,
                        number_color) {
        n = nrow(matrix)
        m = ncol(matrix)

        coord_x = find_coordinates(m, gaps_cols)
        coord_y = find_coordinates(n, gaps_rows)

        x = coord_x$coord - 0.5 * coord_x$size
        y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)

        coord = expand.grid(y = y, x = x)

        res = gList()

        res[["rect"]] =
            rectGrob(x = coord$x,
                     y = coord$y,
                     width = coord_x$size,
                     height = coord_y$size,
                     gp = gpar(fill = matrix, col = border_color))

        if (attr(fmat, "draw")) {
            res[["text"]] =
                textGrob(x = coord$x,
                         y = coord$y,
                         label = fmat,
                         gp = gpar(col = number_color, fontsize = fontsize_number))
        }

        res = gTree(children = res)

        return(res)
    }


draw_colnames <- function(coln, gaps, ...) {
    coord = find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size

    res =
        textGrob(coln,
                 x = x,
                 y = unit(1, "npc") - unit(3, "bigpts"),
                 vjust = 0.5,
                 hjust = 0,
                 rot = 270,
                 gp = gpar(...))

    return(res)
}


draw_rownames <- function(rown, gaps, ...) {
    coord = find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)

    res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))

    return(res)
}


draw_legend <- function(color, breaks, txt.stats, legend.cex, legend, ...) {
    height = min(unit(1 * legend.cex, "npc"), unit(150 * legend.cex, "bigpts")) 

    ## print('*****')
    ## print(height)
    ## print('*****')
    
    ## print(breaks)  

    ## print('dar_legend')  
    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    legend_pos = height * legend_pos + (unit(1, "npc") - height)

    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks = height * breaks + (unit(1, "npc") - height)

    ## print(breaks)  
    ## print(legend_pos[length(legend_pos)])
    ## print(breaks[-length(breaks)])  

    h = breaks[-1] - breaks[-length(breaks)]

    ## print('***** br')
    ## print(breaks)
    ## print('*****')

    ## print(breaks[-length(breaks)]+unit(1, "npc"))

    ## print('***** h')
    ## print(h)
    ## print('*****')
    ## h = cellheight

    rect =
        rectGrob(x = 0,
                 y = breaks[-length(breaks)],
                 width = unit(10, "bigpts"),
                 height = h,
                 hjust = 0,
                 vjust = 0,
                 gp = gpar(fill = color, col = "#FFFFFF00"))
    text =
        textGrob(names(legend),
                 x = unit(14, "bigpts"),
                 y = legend_pos,
                 hjust = 0,
                 gp = gpar(...))

    ## y.stats = breaks[-length(breaks)]-unit(.1, "npc")
    ## stats = rectGrob(x = 0, y = y.stats, width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill ='black'))

    if (!is.na(txt.stats)) {
        ##    rect = rectGrob(x = 0, 
        ##                    y = breaks[-length(breaks)], 
        ##                    width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "red"))

        crlf = strsplit(txt.stats, split = "\\n")[[1]]
        h = length(crlf) / 6

        stats =
            textGrob(txt.stats,
                     x = unit(2, "bigpts"),
                     y = legend_pos[1] - unit(2, "cm"),
                     hjust = 0,
                     gp = gpar(fontface='bold'))

        res = grobTree(rect, text, stats)
    } else
        res = grobTree(rect, text)
    
    return(res)
}


convert_annotations <- function(annotation, annotation_colors) {
    new = annotation
    for (i in 1:ncol(annotation)) {
        a = annotation[, i]
        b = annotation_colors[[colnames(annotation)[i]]]
        if (is.character(a) | is.factor(a)) {
            a = as.character(a)
            if (length(setdiff(a, names(b))) > 0) {
                stop(sprintf("Factor levels on variable %s do not match with annotation_colors",
                             colnames(annotation)[i]))
            }
            new[, i] = b[a]
        }
        else{
            a = cut(a, breaks = 100)
            new[, i] = colorRampPalette(b)(100)[a]
        }
    }
    return(as.matrix(new))
}


draw_annotations <- function(converted_annotations,
                             border_color,
                             gaps, fontsize,
                             horizontal) {
    n = ncol(converted_annotations)
    m = nrow(converted_annotations)

    coord_x = find_coordinates(m, gaps)

    x = coord_x$coord - 0.5 * coord_x$size

    ## y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")

    if (horizontal) {
        coord = expand.grid(x = x, y = y)
        res =
            rectGrob(x = coord$x,
                     y = coord$y,
                     width = coord_x$size,
                     height = unit(fontsize, "bigpts"),
                     gp = gpar(fill = converted_annotations, col = border_color))
    } else {
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a

        coord = expand.grid(y = y, x = x)
        res =
            rectGrob(x = coord$x,
                     y = coord$y,
                     width = unit(fontsize, "bigpts"),
                     height = coord_x$size,
                     gp = gpar(fill = converted_annotations, col = border_color))
    }
    return(res)
}

draw_annotation_names = function(annotations, fontsize, horizontal) {
    n = ncol(annotations)

    x = unit(3, "bigpts")

    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")

    if (horizontal) {
        res = textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
    } else {
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a

        res = textGrob(colnames(annotations),
            x = x,
            y = y,
            vjust = 0.5,
            hjust = 0,
            rot = 270,
            gp = gpar(fontsize = fontsize, fontface = 2))
    }

    return(res)
}


draw_annotation_legend <- function(annotation,
                                   annotation_colors,
                                   border_color,
                                   ...) {
    y = unit(1, "npc")
    text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))

    res = gList()
    for (i in names(annotation)) {
        res[[i]] =
            textGrob(i,
                     x = 0,
                     y = y,
                     vjust = 1,
                     hjust = 0,
                     gp = gpar(fontface = "bold", ...))

        y = y - 1.5 * text_height
        if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
            n = length(annotation_colors[[i]])
            yy = y - (1:n - 1) * 2 * text_height

            res[[paste(i, "r")]] =
                rectGrob(x = unit(0, "npc"),
                         y = yy,
                         hjust = 0,
                         vjust = 1,
                         height = 2 * text_height,
                         width = 2 * text_height,
                         gp = gpar(col = border_color, fill = annotation_colors[[i]]))
            res[[paste(i, "t")]] =
                textGrob(names(annotation_colors[[i]]),
                         x = text_height * 2.4,
                         y = yy - text_height,
                         hjust = 0,
                         vjust = 0.5,
                         gp = gpar(...))

            y = y - n * 2 * text_height

        } else {
            yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
            h = 8 * text_height * 0.25

            res[[paste(i, "r")]] =
                rectGrob(x = unit(0, "npc"),
                         y = yy,
                         hjust = 0,
                         vjust = 1,
                         height = h,
                         width = 2 * text_height,
                         gp =
                             gpar(col = NA,
                                  fill = colorRampPalette(annotation_colors[[i]])(4)))
            
            res[[paste(i, "r2")]] =
                rectGrob(x = unit(0, "npc"),
                         y = y,
                         hjust = 0,
                         vjust = 1,
                         height = 8 * text_height,
                         width = 2 * text_height,
                         gp = gpar(col = border_color))

            txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
            yy = y - c(1, 7) * text_height
            res[[paste(i, "t")]] =
                textGrob(txt,
                         x = text_height * 2.4,
                         y = yy,
                         hjust = 0,
                         vjust = 0.5,
                         gp = gpar(...))
            y = y - 8 * text_height
        }
        y = y - 1.5 * text_height
    }

    res = gTree(children = res)

    return(res)
}


draw_main <- function(text, ...) {
    res = textGrob(text, gp = gpar(fontface = "bold", ...))

    return(res)
}


vplayout <- function(x, y) {
    return(viewport(layout.pos.row = x, layout.pos.col = y))
}


heatmap_motor <- function(matrix,
                          border_color,
                          cellwidth,
                          cellheight,
                          tree_col,
                          tree_row,
                          treeheight_col,
                          treeheight_row,
                          filename,
                          width,
                          height,
                          breaks,
                          color,
                          legend,
                          annotation_row,
                          annotation_col,
                          annotation_colors,
                          annotation_legend,
                          main,
                          fontsize,
                          fontsize_row,
                          fontsize_col,
                          fmat,
                          fontsize_number,
                          number_color,
                          gaps_col,
                          gaps_row,
                          labels_row,
                          labels_col,
                          legend.cex,
                          txt.stats, ...) {
    ## Set layout
    lo =
        lo(coln = labels_col,
           rown = labels_row,
           nrow = nrow(matrix),
           ncol = ncol(matrix),
           cellwidth = cellwidth,
           cellheight = cellheight,
           treeheight_col = treeheight_col,
           treeheight_row = treeheight_row,
           legend = legend,
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           annotation_legend = annotation_legend,
           main = main,
           fontsize = fontsize,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           ...)

    res = lo$gt
    mindim = lo$mindim

    if (!is.na(filename)) {
        if (is.na(height)) {
            height = convertHeight(gtable_height(res), "inches", valueOnly = TRUE)
        }
        if (is.na(width)) {
            width = convertWidth(gtable_width(res), "inches", valueOnly = TRUE)
        }

        ## Get file type.
        
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if (r == -1) stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))

        f =
            switch(ending,
                   pdf = function(x, ...) pdf(x, ...),
                   png = function(x, ...) png(x, units = "in", res = 300, ...),
                   jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
                   jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
                   tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
                   bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
                   stop("File type should be: pdf, png, bmp, jpg, tiff")
                   )

        ## print(sprintf("height:%f width:%f", height, width))

        ## gt = heatmap_motor(matrix, cellwidth = cellwidth,
        ## cellheight = cellheight, border_color = border_color,
        ## tree_col = tree_col, tree_row = tree_row, treeheight_col =
        ## treeheight_col, treeheight_row = treeheight_row, breaks =
        ## breaks, color = color, legend = legend, annotation_col =
        ## annotation_col, annotation_row = annotation_row,
        ## annotation_colors = annotation_colors, annotation_legend =
        ## annotation_legend, filename = NA, main = main, fontsize =
        ## fontsize, fontsize_row = fontsize_row, fontsize_col =
        ## fontsize_col, fmat = fmat, fontsize_number =
        ## fontsize_number, number_color = number_color, labels_row =
        ## labels_row, labels_col = labels_col, gaps_col = gaps_col,
        ## gaps_row = gaps_row, ...)

        f(filename, height = height, width = width)
        gt =
            heatmap_motor(matrix,
                          cellwidth = cellwidth,
                          cellheight = cellheight,
                          border_color = border_color,
                          tree_col = tree_col,
                          tree_row = tree_row,
                          treeheight_col = treeheight_col,
                          treeheight_row = treeheight_row,
                          breaks = breaks,
                          color = color,
                          legend = legend,
                          annotation_col = annotation_col,
                          annotation_row = annotation_row,
                          annotation_colors = annotation_colors,
                          annotation_legend = annotation_legend,
                          filename = NA,
                          main = main,
                          fontsize = fontsize,
                          fontsize_row = fontsize_row,
                          fontsize_col = fontsize_col,
                          fmat = fmat,
                          fontsize_number = fontsize_number,
                          number_color = number_color,
                          labels_row = labels_row,
                          labels_col = labels_col,
                          gaps_col = gaps_col,
                          gaps_row = gaps_row,
                          legend.cex = legend.cex,
                          txt.stats = txt.stats,
                          ...)
        grid.draw(gt)
        dev.off()

        return(gt)
    }

    ## Omit border color if cell size is too small.
    
    if (mindim < 3) border_color = NA

    ## Draw title.
    
    if (!is.na(main)) {
        elem = draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main")
    }

    ## Draw tree for the columns.
    
    if (!is.na(tree_col[[1]][1]) & treeheight_col != 0) {
        elem = draw_dendrogram(tree_col, gaps_col, horizontal = TRUE)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }

    ## Draw tree for the rows.
    
    if (!is.na(tree_row[[1]][1]) & treeheight_row != 0) {
        elem = draw_dendrogram(tree_row, gaps_row, horizontal = FALSE)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }

    ## Draw matrix.
    
    elem = draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "matrix")


    ## Draw colnames.
    
    if (length(labels_col) != 0) {
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
        elem = do.call(draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "col_names")
    }

    ## Giulio
    ##res = gtable_add_grob(
    ##  res, elem, t = 6, l = 3, clip = "off", name = "giulio"
    ##  )
    
    ## Draw rownames.

    if (length(labels_row) != 0) {
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
        elem = do.call(draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "row_names")
    }

    ## Draw annotation tracks on cols.
    
    if (!is.na(annotation_col[[1]][1])) {
        ## Draw tracks.
        
        converted_annotation = convert_annotations(annotation_col, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = TRUE)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", name = "col_annotation")

        ## Draw names.
        
        elem = draw_annotation_names(annotation_col, fontsize, horizontal = TRUE)
        res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", name = "row_annotation_names")
    }

    ## Draw annotation tracks on rows.
    
    if (!is.na(annotation_row[[1]][1])) {
        ## Draw tracks.
        
        converted_annotation = convert_annotations(annotation_row, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = FALSE)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", name = "row_annotation")

        ## Draw names.
        
        elem = draw_annotation_names(annotation_row, fontsize, horizontal = FALSE)
        res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", name = "row_annotation_names")
    }

    ## Draw annotation legend.
    
    annotation = c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na(x[1])))]

    if (length(annotation) > 0 & annotation_legend) {
        elem = draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)

        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
    }

    ## Draw legend.
    
    if (!is.na(legend[1])) {
        elem = draw_legend(color, breaks, txt.stats, legend, legend.cex = legend.cex, fontsize = fontsize, ...)

        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b= 5, clip = "off", name = "legend")
    }

    return(res)
}


generate_breaks <- function(x, n, center = FALSE) {
    if (center) {
        m = max(abs(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))))
        res = seq(-m, m, length.out = n + 1)
    }
    else{
        res = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1)
    }

    return(res)
}

    scale_vec_colours = function(x, col = rainbow(10), breaks = NA) {
        return(col[as.numeric(cut(x, breaks = breaks, include.lowest = TRUE))])
    }

    scale_colours = function(mat, col = rainbow(10), breaks = NA) {
        mat = as.matrix(mat)
        return(matrix(scale_vec_colours(as.vector(mat),
                                        col = col,
                                        breaks = breaks),
                      nrow(mat),
                      ncol(mat),
                      dimnames =
                          list(rownames(mat),
                               colnames(mat))))
    }


cluster_mat <- function(mat, distance, method) {
    if (!(method %in% c("ward.D2",
                        "ward",
                        "single",
                        "complete",
                        "average",
                        "mcquitty",
                        "median",
                        "centroid"))) {
        stop("clustering method must be one of: 'ward', ",
             "'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if (!(distance[1] %in% c("correlation",
                             "euclidean",
                             "maximum",
                             "manhattan",
                             "canberra",
                             "binary",
                             "minkowski"))
        & is(distance, "dist")) {
        stop("distance has to be a dissimilarity structure as produced by dist ",
             "or one measure from the list: ",
             "'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if (distance[1] == "correlation") {
        d = as.dist(1 - cor(t(mat)))
    } else {
        if (is(distance, "dist")) {
            d = distance
        } else {
            d = dist(mat, method = distance)
        }
    }

    return(hclust(d, method = method))
}


scale_rows <- function(x) {
    m = apply(x, 1, mean, na.rm = TRUE)
    s = apply(x, 1, sd, na.rm = TRUE)
    return((x - m) / s)
}


scale_mat <- function(mat, scale) {
    if (!(scale %in% c("none", "row", "column"))) {
        stop("scale argument shoud take values: 'none', 'row' or 'column'")
    }
    mat = switch(scale,
        none = mat,
        row = scale_rows(mat),
        column = t(scale_rows(t(mat))))
    return(mat)
}


generate_annotation_colours <- function(annotation, annotation_colors, drop) {
    if (is.na(annotation_colors)[[1]][1]) {
        annotation_colors = list()
    }
    count = 0
    for (i in 1:length(annotation)) {
        if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
            if (is.factor(annotation[[i]]) & !drop) {
                count = count + length(levels(annotation[[i]]))
            } else {
                count = count + length(unique(annotation[[i]]))
            }
        }
    }

    factor_colors = dscale(factor(1:count), hue_pal(l = 75))

    set.seed(3453)

    cont_counter = 2
    for (i in 1:length(annotation)) {
        if (!(names(annotation)[i] %in% names(annotation_colors))) {
            if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
                n = length(unique(annotation[[i]]))
                if (is.factor(annotation[[i]]) & !drop) {
                    n = length(levels(annotation[[i]]))
                }
                ind = sample(1:length(factor_colors), n)
                annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
                l = levels(as.factor(annotation[[i]]))
                l = l[l %in% unique(annotation[[i]])]
                if (is.factor(annotation[[i]]) & !drop) {
                    l = levels(annotation[[i]])
                }

                names(annotation_colors[[names(annotation)[i]]]) = l
                factor_colors = factor_colors[-ind]
            } else {
                annotation_colors[[names(annotation)[i]]] = brewer_pal("seq", cont_counter)(5)[1:4]
                cont_counter = cont_counter + 1
            }
        }
    }
    return(annotation_colors)
}


kmeans_pheatmap <- function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...) {
    ## Filter data.
    
    if (!is.na(sd_limit)) {
        s = apply(mat, 1, sd)
        mat = mat[s > sd_limit, ]    
    }

    ## Cluster data.
    
    set.seed(1245678)
    km = kmeans(mat, k, iter.max = 100)
    mat2 = km$centers

    ## Compose rownames.
    
    t = table(km$cluster)
    rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)

    ## Draw heatmap.
    
    pheatmap(mat2, ...)
}


find_gaps <- function(tree, cutree_n) {
    v = cutree(tree, cutree_n)[tree$order]
    gaps = which((v[-1] - v[-length(v)]) != 0)
}


#' A function to draw clustered heatmaps.
#' 
#' A function to draw clustered heatmaps where one has better control over some graphical 
#' parameters such as cell size, etc. 
#' 
#' The function also allows to aggregate the rows using kmeans clustering. This is 
#' advisable if number of rows is so big that R cannot handle their hierarchical 
#' clustering anymore, roughly more than 1000. Instead of showing all the rows 
#' separately one can cluster the rows in advance and show only the cluster centers. 
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' This is a modified version of the original pheatmap 
#' (https://cran.r-project.org/web/packages/pheatmap/index.html)
#' edited in accordance with GPL-2.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the 
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors, to certain values. If value is NA then the 
#' breaks are calculated automatically.
#' @param border_color color of cell borders on heatmap, use NA if no border should be 
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values 
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA, 
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in 
#' either the row direction or the column direction, or none. Corresponding values are 
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered,
#' @param cluster_cols boolean values determining if columns should be clustered.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible 
#' values are \code{"correlation"} for Pearson correlation and all the distances 
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible 
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as 
#' \code{\link{hclust}}.
#' @param cutree_rows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the 
#' argument is ignored
#' @param cutree_cols similar to \code{cutree_rows}, but for columns
#' @param treeheight_row the height of a tree for rows, if these are clustered. 
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. 
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation_row data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The 
#' rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotation_col similar to annotation_row, but for columns. 
#' @param annotation deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors list for specifying annotation_row and 
#' annotation_col track colors manually. It is  possible to define the colors 
#' for only some of the features. Check examples for  details.
#' @param annotation_legend boolean value showing if the legend for annotation 
#' tracks should be drawn. 
#' @param drop_levels logical to determine if unused levels are also shown in 
#' the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot 
#' @param legend.cex Default 0.5; determines legend size if \code{legend = TRUE}
#' @param fontsize_row fontsize for rownames (Default: fontsize) 
#' @param fontsize_col fontsize for colnames (Default: fontsize) 
#' @param display_numbers logical determining if the numeric values are also printed to 
#' the cells. If this is a matrix (with same dimensions as original matrix), the contents
#' of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in cells. 
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential 
#' notation (see more in \code{\link{sprintf}}).
#' @param number_color color of the text    
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param gaps_row vector of row indices that show shere to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutree_row}
#'  to see how to introduce gaps to clustered rows. 
#' @param gaps_col similar to gaps_row, but for columns.
#' @param labels_row custom labels for rows that are used instead of rownames.
#' @param labels_col similar to labels_row, but for columns.
#' @param filename file path where to save the picture. Filetype is decided by 
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param txt.stats By default, shows a summary statistics for shown data (n,m, |G| and |P|)
#' @param \dots graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @return 
#' Invisibly a list of components 
#' \itemize{
#'     \item \code{tree_row} the clustering of rows as \code{\link{hclust}} object 
#'     \item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was 
#' specified 
#' }
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Create test matrix
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' 
#' # Draw heatmaps
#' pheatmap(test)

#' @export pheatmap
#' @importFrom grid unit textGrob gpar unit.c viewport convertWidth
#' @importFrom grid convertHeight gList gTree rectGrob grobTree polylineGrob
#' @importFrom grid grid.draw grid.pretty grid.newpage
#' @importFrom scales dscale hue_pal brewer_pal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats as.dist cor hclust dist cutree sd kmeans sd
#' @importFrom grDevices colorRampPalette pdf png jpeg tiff bmp dev.off rainbow
#' @importFrom graphics strwidth
#' 
pheatmap <- function(mat,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                     kmeans_k = NA,
                     breaks = NA,
                     border_color = "grey60",
                     cellwidth = NA,
                     cellheight = NA,
                     scale = "none",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     cutree_rows = NA,
                     cutree_cols = NA,
                     treeheight_row = ifelse(cluster_rows, 50, 0),
                     treeheight_col = ifelse(cluster_cols, 50, 0),
                     legend = TRUE,
                     legend_breaks = NA,
                     legend_labels = NA,
                     annotation_row = NA,
                     annotation_col = NA,
                     annotation = NA,
                     annotation_colors = NA,
                     annotation_legend = TRUE,
                     drop_levels = TRUE,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     main = NA,
                     fontsize = 10,
                     fontsize_row = fontsize,
                     fontsize_col = fontsize,
                     display_numbers = FALSE,
                     number_format = "%.2f",
                     number_color = "grey30",
                     fontsize_number = 0.8 * fontsize,
                     gaps_row = NULL,
                     gaps_col = NULL,
                     labels_row = NULL,
                     labels_col = NULL,
                     filename = NA,
                     width = NA,
                     height = NA,
                     silent = FALSE,
                     legend.cex = 1,
                     txt.stats = NA,
                     ...) {
    ## Set labels.
    
    if (is.null(labels_row)) {
        labels_row = rownames(mat)
    }
    if (is.null(labels_col)) {
        labels_col = colnames(mat)
    }

    ## Preprocess matrix.
    
    mat = as.matrix(mat)
    if (scale != "none") {
        mat = scale_mat(mat, scale)
        if (is.na(breaks)) {
            breaks = generate_breaks(mat, length(color), center = TRUE)
        }
    }

    ## Kmeans.
    
    if (!is.na(kmeans_k)) {
        ## Cluster data.
        
        km = kmeans(mat, kmeans_k, iter.max = 100)
        mat = km$centers

        ## Compose rownames.
        
        t = table(km$cluster)
        labels_row = sprintf("Cluster: %s Size: %d", names(t), t)
    } else {
        km = NA
    }

    ## Format numbers to be displayed in cells.
    
    if (is.matrix(display_numbers) | is.data.frame(display_numbers)) {
        if (nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)) {
            stop("If display_numbers provided as matrix, its dimensions have to match with mat")
        }

        display_numbers = as.matrix(display_numbers)
        fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
        fmat_draw = TRUE
    } else {
        if (display_numbers) {
            fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = TRUE
        } else {
            fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = FALSE
        }
    }

    ## Do clustering.
    
    if (cluster_rows) {
        tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
        mat = mat[tree_row$order, , drop = FALSE]
        fmat = fmat[tree_row$order, , drop = FALSE]
        labels_row = labels_row[tree_row$order]
        if (!is.na(cutree_rows)) {
            gaps_row = find_gaps(tree_row, cutree_rows)
        } else {
            gaps_row = NULL
        }
    } else {
        tree_row = NA
        treeheight_row = 0
    }

    if (cluster_cols) {
        tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
        mat = mat[, tree_col$order, drop = FALSE]
        fmat = fmat[, tree_col$order, drop = FALSE]
        labels_col = labels_col[tree_col$order]
        if (!is.na(cutree_cols)) {
            gaps_col = find_gaps(tree_col, cutree_cols)
        } else {
            gaps_col = NULL
        }
    } else {
        tree_col = NA
        treeheight_col = 0
    }

    attr(fmat, "draw") = fmat_draw

    ## Colors and scales.
    
    if (!is.na(legend_breaks[1]) & !is.na(legend_labels[1])) {
        if (length(legend_breaks) != length(legend_labels)) {
            stop("Lengths of legend_breaks and legend_labels must be the same")
        }
    }

    if (is.na(breaks[1])) {
        breaks = generate_breaks(as.vector(mat), length(color))
    }
    
    if (legend & is.na(legend_breaks[1])) {
        legend = grid.pretty(range(as.vector(breaks)))
        names(legend) = legend
    } else if (legend & !is.na(legend_breaks[1])) {
        legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]

        if (!is.na(legend_labels[1])) {
            legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
            names(legend) = legend_labels
        } else {
            names(legend) = legend
        }
    } else {
        legend = NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)

    ## Preparing annotations.
    
    if (is.na(annotation_col[[1]][1]) & !is.na(annotation[[1]][1])) {
        annotation_col = annotation
    }
    
    ## Select only the ones present in the matrix.
    
    if (!is.na(annotation_col[[1]][1])) {
        annotation_col = annotation_col[colnames(mat), , drop = FALSE]
    }

    if (!is.na(annotation_row[[1]][1])) {
        annotation_row = annotation_row[rownames(mat), , drop = FALSE]
    }

    annotation = c(annotation_row, annotation_col)
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na(x[1])))]
    if (length(annotation) != 0) {
        annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
    } else {
        annotation_colors = NA
    }

    if (!show_rownames) {
        labels_row = NULL
    }

    if (!show_colnames) {
        labels_col = NULL
    }

    ## Draw heatmap
    gt =
        heatmap_motor(mat,
                      border_color = border_color,
                      cellwidth = cellwidth,
                      cellheight = cellheight,
                      treeheight_col = treeheight_col,
                      treeheight_row = treeheight_row,
                      tree_col = tree_col,
                      tree_row = tree_row,
                      filename = filename,
                      width = width,
                      height = height,
                      breaks = breaks,
                      color = color,
                      legend = legend,
                      annotation_row = annotation_row,
                      annotation_col = annotation_col,
                      annotation_colors = annotation_colors,
                      annotation_legend = annotation_legend,
                      main = main,
                      fontsize = fontsize,
                      fontsize_row = fontsize_row,
                      fontsize_col = fontsize_col,
                      fmat = fmat,
                      fontsize_number = fontsize_number,
                      number_color = number_color,
                      gaps_row = gaps_row,
                      gaps_col = gaps_col,
                      labels_row = labels_row,
                      labels_col = labels_col,
                      legend.cex = legend.cex,
                      txt.stats = txt.stats,
                      ...)

    if (is.na(filename) & !silent) {
        grid.newpage()
        grid.draw(gt)
    }

    invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt))
}


#' tronco.pattern.plot : plot a genotype
#' 
#' @title tronco.pattern.plot
#' @param x A TRONCO compliant dataset
#' @param group A list of events (see as.events() for details)
#' @param to A target event
#' @param gap.cex cex parameter for gap
#' @param legend.cex cex parameter for legend
#' @param label.cex cex parameter for label
#' @param title title
#' @param mode can be 'circos' or 'barplot'
#' @export tronco.pattern.plot
#' @importFrom circlize circos.clear circos.par chordDiagram circos.text
#' @importFrom circlize circos.trackPlotRegion get.cell.meta.data
#' @importFrom graphics layout par barplot plot.new legend
#' @importFrom grDevices rgb col2rgb
#' 
tronco.pattern.plot <- function(x,
                         group=as.events(x),
                         to,
                         gap.cex = 1.0,
                         legend.cex = 1.0,
                         label.cex = 1.0,
                         title=paste(to[1], to[2]),
                         mode = 'barplot') {
    ## Keys.
    
    events = as.events(x)

    keys = NULL
    for (i in 1:nrow(group))
        keys =
            c(keys,
              rownames(as.events(x,
                                 genes = group[i, 'event'],
                                 types = group[i, 'type'])))

    cat('Group:\n')
    events.names = events[keys, , drop = FALSE]
    print(events.names)

    cat('Group tested against:', to[1], to[2], '\n')

    ## HARD exclusivity: 1 for each pattern element
    ## CO-OCCURRENCE: 1
    ## SOFT exclusivity: 1
    ## OTHERS: 1
    
    matrix = matrix(0, nrow = length(keys) + 3, ncol = 1)
    rownames(matrix) = c(keys, 'soft', 'co-occurrence', 'other')
    ## colnames(matrix) = paste(to, collapse = ':')
    colnames(matrix) = to[1]

    to.samples = which.samples(x, gene=to[1], type=to[2])
    cat('Pattern conditioned to samples:', paste(to.samples, collapse = ', '), '\n')

    pattern.samples = list()
    hard.pattern.samples = list()
    for (i in 1:length(keys)) {
        ## Samples
        samples = which.samples(x, gene= events.names[i, 'event'], type= events.names[i, 'type'])

        ## Pattern samples
        pattern.samples = append(pattern.samples, list(samples))

        ## Other negative samples
        negative.samples = list()
        for (j in 1:length(keys))
            if (i != j)
                negative.samples =
                    append(negative.samples, 
                           list(which.samples(x,
                                              gene = events.names[j, 'event'],
                                              type = events.names[j, 'type'],
                                              neg = TRUE)))

        pattern.negative = Reduce(intersect, negative.samples)
        hard.pattern.samples = append(hard.pattern.samples, list(intersect(samples, pattern.negative)))   
    }
    names(pattern.samples) = keys
    names(hard.pattern.samples) = keys

    ## print(hard.pattern.samples)
    
    ## print('tutti:')
    ## print(unlist(pattern.samples))

    ## CO-OCCURRENCES.
    
    pattern.co.occurrences = Reduce(intersect, pattern.samples)
    co.occurrences = intersect(to.samples, pattern.co.occurrences)
    matrix['co-occurrence', ] = length(co.occurrences)
    cat('Co-occurrence in #samples: ', matrix['co-occurrence', ], '\n')

    ## HARD EXCLUSIVITY.
    
    for (i in 1:length(keys))
        matrix[keys[i], ] = length(intersect(to.samples, hard.pattern.samples[[keys[i]]])) 
    cat('Hard exclusivity in #samples:', matrix[keys, ], '\n')  

    ## OTHER.
    
    union = unique(unlist(pattern.samples))
    union = setdiff(as.samples(x), union)
    matrix['other', ] = length(intersect(to.samples, union))
    cat('Other observations in #samples:', matrix['other', ], '\n') 

    ## SOFT EXCLUSIVITY.
    
    matrix['soft', ] = length(to.samples) - colSums(matrix)
    cat('Soft exclusivity in #samples:', matrix['soft', ], '\n')  

    ## Choose colors according to event type.
    
    sector.color = rep('gray', nrow(matrix) + 1) 
    link.color = rep('gray', nrow(matrix)) 

    names(sector.color) = c(rownames(matrix), colnames(matrix))
    names(link.color) = rownames(matrix)

    add.alpha <- function(col, alpha = 1) {
        if (missing(col)) stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
    }

    ## Link colors - event types .
    
    for (i in 1:length(keys))
        link.color[keys[i]] = as.colors(x)[events.names[keys[i], 'type' ]]

    ## print(link.color)

    link.color['soft'] = 'orange'
    link.color['co-occurrence'] = 'darkgreen'


    idx.max = which(matrix == max(matrix))
    link.style = matrix(0, nrow=nrow(matrix), ncol=ncol(matrix))
    rownames(link.style) = rownames(matrix)
    colnames(link.style) = colnames(matrix)
    link.style[idx.max, 1] = 5
    
    ## print(link.style)

    ## link.color[idx.max] = add.alpha(link.color[idx.max], .3)
    ## print(link.color)


    ## Sector colors - event types.
    
    sector.color[1:length(keys)] = 'red' # XOR
    sector.color['soft'] = 'orange'      # OR
    sector.color['co-occurrence'] = 'darkgreen' # AND
    sector.color[colnames(matrix)] = as.colors(x)[as.events(x, genes = to[1], types=to[2])[, 'type' ]]

    ## print(link.color)
    ## print(sector.color)

    ## Informative labels
    for (i in 1:length(keys)) {
        ## rownames(matrix)[i] = paste(paste(rep(' ', i), collapse = ''), events.names[i, 'event' ])
        if (nevents(x, genes = events.names[i, 'event' ]) > 1)
            rownames(matrix)[i] = paste(paste(rep(' ', i), collapse = ''), events.names[i, 'event' ])
        else rownames(matrix)[i] = events.names[i, 'event' ]

        names(sector.color)[i] = rownames(matrix)[i]    
    }

    if ( mode == 'circos') { 

        cat('Circlize matrix.\n')
        print(matrix)

        circos.clear()

        gaps = c(rep(2, length(keys) - 2), rep(15 * gap.cex, 4), rep(40 * gap.cex, 2))

        circos.par(gap.degree = gaps)

        chordDiagram(matrix, 
                     grid.col = sector.color,
                     annotationTrack = "grid", 
                     preAllocateTracks = list(track.height = 0.3),
                     row.col = link.color,
                     link.border = 'black',
                     link.lty = link.style,
                     link.lwd = 0.3,
                     reduce = 0
                     )

        circos.trackPlotRegion(track.index = 1, 
                               panel.fun = function(x, y) {
                                   xlim = get.cell.meta.data("xlim")
                                   ylim = get.cell.meta.data("ylim")
                                   sector.name = get.cell.meta.data("sector.index")
                                   circos.text(mean(xlim),
                                               cex = 1.0 * label.cex,
                                               ylim[1],
                                               sector.name,
                                               facing = "clockwise",
                                               niceFacing = TRUE,
                                               adj = c(0, 0.5)) 
                               },
                               bg.border = NA)
        
    } else { # barplot
        
        layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), heights = c(4, 1))
        par(mai=rep(0.5, 4))

        basic.unit = 2
        widths =
            c(rep(basic.unit, length(keys)), 
              rep(basic.unit * 2, 3)
              )

        spaces =
            c(rep(0, length(keys)),
              3, 
              rep(.5, 2)
              )

        print(matrix)

        ## barplot(matrix[, 1], widths, space = 0)

        rownames(matrix)[length(keys) + 1] = '2 or more\n(soft-exclusivity)'
        rownames(matrix)[length(keys) + 2] = 'all together\n(co-occurrence)'

        rownames(matrix)[nrow(matrix)] = 'none of the\nevents'

        summary = matrix[ (length(keys) + 1):nrow(matrix), 1, drop = FALSE]
        summary = rbind( sum(matrix[1:length(keys),]), summary)
        rownames(summary) = NULL
        colnames(summary) = paste(to[1], to[2])
        print(summary)  

        ## Plot, but suppress the labels.

        midpts =
            barplot(summary,
                    5,
                    col = c('red', 'orange', 'darkgreen', 'gray'),
                    horiz = FALSE,
                    space = 1,
                    las = 1,
                    main =
                        paste0('Combination of group events \n in ',
                               sum(summary),
                               ' samples with ',
                               to[1],
                               ' ',
                               to[2]),
                    cex.main = .6,
                    cex.names = .5,
                    cex.axis = .5 )

        exclus = matrix[1:length(keys), 1, drop = FALSE]
        ## print(exclus)


        midpts =
            barplot(exclus[, 1],
                    ## widths,
                    col = link.color[1:length(keys)],
                    horiz = TRUE,
                    space = 1,
                    las = 2,
                    main = paste0('Observations supporting \n hard-exclusivity'),
                    cex.main = 0.6,
                    ## xlab = 'number of observations (given KRAS)',
                    cex.names = 0.5,
                    cex.axis = 0.5 )

        par(mai = c(0, 0, 0, 0))
        plot.new()
        
        }

        events.legend.pch = rep(19, length(unique(events.names[, 'type'])))

        legend("topleft",
               cex = 0.6  * legend.cex,
               pt.cex = 1,
               title = expression(bold('Table of observations')),
               horiz=FALSE,
               bty='n',
               legend =
                   c(paste(sum(matrix[1:length(keys) ,]), 'with 1 event (hard exclusivity)'),
                     paste(matrix[length(keys) + 1, ],  'with 2 or more events'),
                     paste(matrix[length(keys) + 2, ], 'with all events (co-occurrence)'),
                     paste(matrix[length(keys) + 3, ], 'with no events')
                     ),
               fill = c('red', 'orange', 'darkgreen', 'gray')
               )

        group.legend = apply(group, 1, paste, collapse = ' in ')

        legend('top',
               cex = .6 * legend.cex,
               pt.cex = 1,
               title = expression(bold('Input group')),
               bty = 'n',
               legend = group.legend
               )

        legend(x = 'topright',
               legend = unique(events.names[, 'type']),
               title = expression(bold('Events type')),
               bty = 'n',
               inset  = +.05,
               cex = .6 * legend.cex,
               pt.cex = 1,
               pch = events.legend.pch,
               col = as.colors(x)[unique(events.names[, 'type'])]
               ## pt.bg = pt_bg
               )
    
}

#### end of file -- visualization.R
