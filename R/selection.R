#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2017, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' select a subset of the input genotypes 'x'. Selection can be done 
#' by frequency and gene symbols.
#' @title events.selection
#'
#' @examples
#' data(test_dataset)
#' dataset = events.selection(test_dataset, 0.3)
#'
#' @param x A TRONCO compliant dataset.
#' @param filter.freq [0,1] value which constriants the minimum frequence of selected events
#' @param filter.in.names gene symbols which will be included
#' @param filter.out.names gene symbols which will NOT be included
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant dataset.
#' @export events.selection
#' 
events.selection <- function(x,
                             filter.freq = NA,
                             filter.in.names = NA,
                             filter.out.names = NA,
                             silent = FALSE) {

    is.compliant(x, err.fun='events.selection: input')
    dataset = x$genotypes

    if (!silent) {
        cat('*** Events selection: #events = ',
            nevents(x),
            ', #types = ',
            ntypes(x))

        cat(' Filters freq|in|out = {', 
            !is.na(filter.freq), ', ',
            !any(is.na(filter.in.names)), ', ',
            !any(is.na(filter.out.names)), '}')
    }


    if (is.na(filter.out.names) && is.na(filter.in.names) && is.na(filter.freq)) {
        return(x)
    }

    valid = rep(FALSE, ncol(x$genotypes))

    if (!is.na(filter.freq)) {

        if (!silent) {
            cat('\nMinimum event frequency: ',
                filter.freq,
                ' (',
                round(nsamples(x) * filter.freq, 0),
                ' alterations out of ',
                nsamples(x),
                ' samples).\n')
        }
        x = enforce.numeric(x)


        #flush.console()
        #pb = txtProgressBar(1, nevents(x), style = 3)
        
        for (i in 1:nevents(x)) {   
            #setTxtProgressBar(pb, i)
            mut.freq = sum(x$genotypes[,i])/nsamples(x)
            valid[i] = mut.freq > filter.freq
            if (!silent) {
                cat('.')
            }
        }
        
        #close(pb)

        if (!silent) {
            cat('\nSelected ',
                nrow(as.events(x)[valid, ]),
                ' events.\n')
        }
    }

    if (!any(is.na(filter.in.names))) {
        shown = min(5, length(filter.in.names))

        if (!silent) {
            cat('\n[filter.in] Genes hold: ', 
                paste(filter.in.names[1:shown], collapse=', '),
                ' ... ')
        }

        colnames = which(x$annotations[,2] %in% filter.in.names,
            arr.ind = TRUE)

        k = unique(x$annotations[which(x$annotations[ ,'event'] %in% filter.in.names, 
                                    arr.ind = TRUE),
                                    'event'])

        if (!silent) {
            cat(' [', length(k), '/', length(filter.in.names), ' found].')
        }

        valid[colnames] = TRUE
    }

    if (!any(is.na(filter.out.names))) {
        shown = min(5, length(filter.out.names))

        if (!silent) {
            cat('\n[filter.out] Genes dropped: ', 
                paste(filter.out.names[1:shown], collapse=', '),
                ' ... ')
        }

        colnames = which(x$annotations[,2] %in% filter.out.names,
            arr.ind = TRUE)

        if (!silent) {
            cat(' [',
                length(colnames),
                '/',
                length(filter.out.names),
                ' found].')
        }
        valid[colnames] = FALSE
    }

    y = list()
    y$genotypes = x$genotypes[, valid, drop = FALSE]  

    y$annotations = as.matrix(x$annotations[valid, , drop = FALSE])
    colnames(y$annotations) = c('type', 'event')
    rownames(y$annotations) = colnames(y$genotypes)

    y$types = as.matrix(x$types[unique(y$annotations[,1]), 1])
    colnames(y$types) = c('color')
    rownames(y$types) = unique(y$annotations[,1])

    if (!is.null(x$stages)) y$stages=x$stages
    is.compliant(x, err.fun='events.selection: output')

    if (!silent) {
        cat('\nSelected ', nevents(y), ' events, returning.\n')
    }

    return(y)
}


#' Return the first n recurrent events
#' @title rank.recurrents
#'
#' @examples
#' data(test_dataset)
#' dataset = rank.recurrents(test_dataset, 10)
#'
#' @param x A TRONCO compliant dataset.
#' @param n The number of events to rank 
#' @return the first n recurrent events
#' @export rank.recurrents
#' 
rank.recurrents <- function(x, n) {
    is.compliant(x)
    x = enforce.numeric(x)    

    if(n <= 0) {
        stop('Rank value (n) should be positive.')
    }

    ## Sum columns
    sums = colSums(x$genotypes)

    ## Get the names of the first n ranked
    sorted = sort(sums, decreasing = TRUE)

    scores = unique(sorted)
    
    l = length(scores)
    if(n >l) {
        warning(paste0('Rank contains ',
                       l,
                       ' unique entries, using n = ',
                       l,
                       ' instead of n = ',
                       n))
    }

    n = min(n, length(scores))
    scores = scores[1:n]

    sorted = sorted[which(sorted >= min(scores))]

    max = names(sorted[which(sorted == max(scores))])
    min = names(sorted[which(sorted == min(scores))])

    cat(paste0('Most recurrent(s): ',
               paste(as.events(x)[max, 'event'], collapse=', '),
               ' (',
               (max(scores)),
               ' hits).\n' ))
    cat(paste0(n, '-th recurrent(s): ',
               paste(as.events(x)[min, 'event'], collapse=', '),
               ' (',
               (min(scores)),
               ' hits).\n' ))

    order = names(sorted)
    genes = as.events(x)[order, 'event']

    return(as.vector(genes))
}

#' Filter a dataset based on selected samples id
#' @title samples.selection
#'
#' @examples
#' data(test_dataset)
#' dataset = samples.selection(test_dataset, c('patient 1', 'patient 2'))
#'
#' @param x A TRONCO compliant dataset.
#' @param samples A list of samples
#' @return A TRONCO compliant dataset.
#' @export samples.selection
#' 
samples.selection <- function(x, samples) {
    is.compliant(x, 'Input:')

    missing = setdiff(samples, as.samples(x))
    if (length(missing) > 0) {
        warning(paste('Missing samples: ', paste(missing, collapse=', ')))
    }

    delete = setdiff(as.samples(x), samples)
    return(delete.samples(x, delete))
}

#### end of file -- selection.R
