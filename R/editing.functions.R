#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# Internal function

#' Verify if the input data are consolidate, i.e., if there are events
#' with 0 or 1 probability or indistinguishable in terms of
#' observations
#' 
#' @title consolidate.data
#' 
#' @examples
#' data(test_dataset)
#' consolidate.data(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @param print A boolean value stating whether to print of not the summary
#' @return The list of any 0 probability, 1 probability and indistinguishable.
#' @export consolidate.data
#' 
consolidate.data <- function(x, print = FALSE) {

    is.compliant(x)

    ind = list()
    zeros = list()
    ones = list()

    ## Get the list of indistinguishible events for each entry.

    duplicated.names =
        unlist(lapply(colnames(x$genotypes),
                      function(gene)
                          paste(as.vector(which(apply(t(x$genotypes),
                                                      1,
                                                      function(y)
                                                          all(y == as.vector(x$genotypes[ , gene]))))),
                                collapse = "-")))

    ## Get the unique set of indistinguishible events.

    duplicated.events =
        unique(duplicated.names[which(unlist(lapply(duplicated.names,
                                                    ## The next
                                                    ## deserves to be
                                                    ## sent to dailywtf!
                                                    function(x)
                                                        if (length(unlist(strsplit(x, "-"))) > 1)
                                                            TRUE
                                                        else
                                                            FALSE)))])

    ## Evaluate the indistinguishible events.
    if (length(duplicated.events) > 0) {
        for (i in duplicated.events) {
            ev = as.numeric(unlist(strsplit(i, "-")))

            ind.pool = as.events(x)[ev,]

            if (any(x$genotypes[, ev[1]] != 0) && any(x$genotypes[, ev[1]] != 1)) {
                ind = append(ind, list(ind.pool))
            }

            if (print) {
                cat('\nIndistinguishable events:\n')

                print(ind.pool)
                cat('Total number of events for these genes: ')
                cat(paste(nevents(x, as.events(x, genes=c(ind.pool)))), '\n')
            }
        }
    }

    ## Evaluate any event at 0 or 1 probability.

    for (i in 1:ncol(x$genotypes)) {

        ev.pool = as.events(x)[i,]

        if (sum(x$genotypes[, i]) == 0) {
            zeros = append(zeros, list(ev.pool))
            if (print) {

                cat('\nIndistinguishable events:\n')

                print(ev.pool)
                cat('Total number of events for these genes: ')
                cat(paste(nevents(x, as.events(x, genes = c(ev.pool)))), '\n')

            }
        }
        if (sum(x$genotypes[, i]) == nrow(x$genotypes)) {
            ones = append(ones, list(ev.pool))
            if (print) {
                cat('\nEvents with no alterations across samples:\n')

                print(ev.pool)
                cat('Total number of events for these genes: ')
                cat(paste(nevents(x, as.events(x, genes=c(ev.pool)))), '\n')
            }
        }
    }

    ret = NULL
    ret$indistinguishable = ind
    ret$zeroes = zeros
    ret$ones = ones

    return(ret)
}


#' Annotate a description on the selected dataset
#' @title annotate.description
#'
#' @examples
#' data(test_dataset)
#' annotate.description(test_dataset, 'new description')
#'
#' @param x A TRONCO compliant dataset.
#' @param label A string
#' @return A TRONCO compliant dataset.
#' @export annotate.description
#' 
annotate.description <- function(x, label) {
    if (as.description(x) != "")
        warning(paste('Old description substituted: ', as.description(x), '.'))

    x$name = label
    return(x)
}


#' Annotate stage information on the selected dataset
#' @title annotate.stages
#'
#' @examples
#' data(test_dataset)
#' data(stage)
#' test_dataset = annotate.stages(test_dataset, stage)
#' as.stages(test_dataset)
#' 
#' @param x A TRONCO compliant dataset.
#' @param stages A list of stages. Rownames must match samples list of x
#' @param match.TCGA.patients Match using TCGA notations (only first 12 characters)
#' @return A TRONCO compliant dataset.
#' @export annotate.stages
#' 
annotate.stages <- function(x, stages, match.TCGA.patients = FALSE) {
    if (is.null(rownames(stages))) {
        stop('Stages have no rownames - will not add annotation.')
    }

    samples = as.samples(x)

    ## Just for temporary - will be shortened to make a simple check...
    
    samples.temporary = samples
    if (match.TCGA.patients) {
        samples.temporary = substring(samples.temporary, 0, 12)
    }

    if (!any(samples.temporary %in% rownames(stages)))
        stop('There are no stages for samples in input dataset - will not add annotation.')

    ## Notify if something gets lost.
    
    if (has.stages(x)) {
        warning('Stages in input dataset overwritten.')
    }

    ## Actual stages.
    
    x$stages = data.frame(row.names = samples, stringsAsFactors = FALSE)
    x$stages[ , 'stage'] = as.character(NA)

    for (i in 1:nsamples(x)) {
        if (!match.TCGA.patients) {
            x$stages[i, ] = as.character(stages[as.samples(x)[i], ])
        } else {
            
            ## Potential match if x's samples are long TCGA barcodes
            ## and stages are TCGA patients barcdodes (short).
            
            short.name = substring(as.samples(x)[i], 0 , 12)
            x$stages[i, 'stage'] = as.character(stages[short.name, ])
        }
    }

    count.na = is.na(x$stages)
    if (any(count.na)) {
        warning(paste(length(which(count.na)),
                      ' missing stages were added as NA.'))
    }

    return(x)
}


#' Change the color of an event type
#' @title change.color
#'
#' @examples
#' data(test_dataset)
#' dataset = change.color(test_dataset, 'ins_del', 'red')
#'
#' @param x A TRONCO compliant dataset.
#' @param type An event type
#' @param new.color The new color (either HEX or R Color)
#' @return A TRONCO complian dataset.
#' @export change.color
#' 
change.color <- function(x, type, new.color) {
    is.compliant(x)
    if (type %in% as.types(x)) {
        x$types[type, ] = new.color
    } else {
        stop('type: \"', type, '\" not in dataset')
    }

    is.compliant(x)
    return(x)
}


#' Rename an event type
#' @title rename.type
#'
#' @examples
#' data(test_dataset)
#' test_dataset = rename.type(test_dataset, 'ins_del', 'deletion')
#'
#' @param x A TRONCO compliant dataset.
#' @param old.name The type of event to rename.
#' @param new.name The new name
#' @return A TRONCO complian dataset.
#' @export rename.type
#' 
rename.type <- function(x, old.name, new.name) {
    is.compliant(x, 'rename.type: input dataset')
    types = as.types(x)

    if (old.name == 'Pattern'
        || new.name == 'Pattern') {
        stop("Pattern is a reserved keyword in TRONCO")
    }

    if (old.name == new.name) {
        return(x)
    }


    if (old.name %in% types) {
        x$annotations[which(x$annotations[ , 'type'] == old.name), 'type'] = new.name
        if (! new.name %in% types) {
            rownames(x$types)[which(rownames(x$types) == old.name)] = new.name
        } else {
            x$types =
                x$types[!rownames(x$types) %in% list(old.name), , drop = FALSE]
        }
    } else {
        stop(paste(old.name, 'not in as.types(x)'))
    }

    is.compliant(x, err.fun = 'rename.type: output')
    return(x)
}


#' Rename a gene
#' @title rename.gene
#'
#' @examples
#' data(test_dataset)
#' test_dataset = rename.gene(test_dataset, 'TET2', 'gene x')
#'
#' @param x A TRONCO compliant dataset.
#' @param old.name The name of the gene to rename.
#' @param new.name The new name
#' @return A TRONCO complian dataset.
#' @export rename.gene
#' 
rename.gene <- function(x, old.name, new.name) {
    ## If is compliant x
    is.compliant(x)

    if (old.name %in% as.genes(x)) {
        x$annotations[ which(x$annotations[,'event'] == old.name), 'event' ] = new.name

    } else {
        stop(paste(old.name, 'not in as.genes(x)'))
    }

    is.compliant(x)
    return(x)
}


#' Delete an event type
#' @title delete.type
#'
#' @examples
#' data(test_dataset)
#' test_dataset = delete.type(test_dataset, 'Pattern')
#'
#' @param x A TRONCO compliant dataset.
#' @param type The name of the type to delete.
#' @return A TRONCO complian dataset.
#' @export delete.type
#' 
delete.type <- function(x, type) {
    ## if is compliant x
    is.compliant(x)

    if (has.model(x)) {
        stop("There's a reconstructed model, a type cannot be deleted now. \nUse delete.model()")
    }

    for (pattern in as.patterns(x)) {
        if (type %in% as.types.in.patterns(x, patterns=pattern)) {
            stop('Found type \"',
                 type,
                 '\" in pattern \"',
                 pattern,
                 '\". Delete that pattern first.\n')
        }
    }

    if (type %in% as.types(x)) {
        events = as.events(x, types = setdiff(as.types(x), type))

        x$genotypes = as.genotypes(x)[ , rownames(events), drop = FALSE]
        x$annotations = events
        x$types = x$types[which(rownames(x$types) != type), , drop = FALSE]
    } else {
        stop(paste(type, 'not in as.types(x)'))
    }

    if (type == 'Pattern') {
        x$hypotheses = NA
    }

    is.compliant(x)
    return(x)
}


#' Delete a gene
#' @title delete.gene
#'
#' @examples
#' data(test_dataset)
#' test_dataset = delete.gene(test_dataset, 'TET2')
#'
#' @param x A TRONCO compliant dataset.
#' @param gene The name of the gene to delete.
#' @return A TRONCO complian dataset.
#' @export delete.gene
#' 
delete.gene <- function(x, gene) {
    ## if is compliant x
    is.compliant(x, 'delete:gene: input')

    if (has.model(x)) {
        stop("There's a reconstructed model, a type cannot be deleted now. \nUse delete.model()")
    }

    if (all(gene %in% as.genes(x))) {
        for (pattern in as.patterns(x)) {
            for (g in gene) {
                if (g %in% as.genes.in.patterns(x, patterns=pattern)) {
                    stop('Found gene \"',
                         g,
                         '\" in pattern \"',
                         pattern,
                         '\". Delete that pattern first.\n')
                }
            }
        }


        drops = rownames(as.events(x, genes = gene))
        x$genotypes = x$genotypes[ , -which(colnames(x$genotypes) %in% drops)]
        x$annotations = x$annotations[ -which(rownames(x$annotations) %in% drops), ]

        ## TO DO: something better than this t(t(...))
        x$types = x$types[which(rownames(x$types) %in% unique(x$annotations[,"type"])), , drop = FALSE]
        colnames(x$types) = 'color'
    } else {
        stop(paste(gene[!(gene %in% as.genes(x))], collapse= ','),
             'are not in the dataset -- as.genes(x).')
    }

    is.compliant(x, 'delete:gene: output')
    return(x)
}


#' Delete an event from the dataset
#' @title delete.event
#'
#' @examples
#' data(test_dataset)
#' test_dataset = delete.event(test_dataset, 'TET2', 'ins_del')
#'
#' @param x A TRONCO compliant dataset.
#' @param gene The name of the gene to delete.
#' @param type The name of the type to delete.
#' @return A TRONCO complian dataset.
#' @export delete.event
#' 
delete.event <- function(x, gene, type) {

    is.compliant(x, 'delete.event: input')

    if (has.model(x)) {
        stop("There's a reconstructed model, a type cannot be deleted now. \nUse delete.model()")
    }

    if (type == 'Pattern') {
        stop("Pattern is a reserved keyword in TRONCO. See delete.pattern")
    }

    for (pattern in as.patterns(x)) {
        events = as.events.in.patterns(x, patterns=pattern)
        if (length(which(events[,'type'] == type & events[,'event'] == gene)) > 0) {
            stop('Found event \"(',
                 gene,
                 ', ',
                 type,
                 ')\" in pattern \"',
                 pattern,
                 '\". Delete that pattern first.\n')
        }
    }


    if (all(c(type, gene) %in% as.events(x))) {
        drops = rownames(as.events(x, genes = gene, types = type))
        x$genotypes = x$genotypes[, -which( colnames(x$genotypes) %in% drops ), drop = FALSE]
        x$annotations = x$annotations[ -which (rownames(x$annotations) %in% drops), , drop = FALSE]

        ## TO DO: something better than this t(t(...))
        x$types = x$types[ which(rownames(x$types) %in% unique(x$annotations[,"type"])), , drop = FALSE]
    } else {
        stop(paste(type, gene, ' - not in as.events(x)'))
    }

    is.compliant(x, 'delete:gene: output')
    return(x)
}


#' Delete an hypothesis from the dataset based on a selected event.
#' Check if the selected event exist in the dataset and delete his associated hypothesis
#' @title delete.hypothesis
#'
#' @examples
#' data(test_dataset)
#' delete.hypothesis(test_dataset, event='TET2')
#' delete.hypothesis(test_dataset, cause='EZH2')
#' delete.hypothesis(test_dataset, event='XOR_EZH2')
#'
#' @param x A TRONCO compliant dataset.
#' @param event Can be an event or pattern name
#' @param cause Can be an event or pattern name
#' @param effect Can be an event or pattern name
#' @return A TRONCO complian dataset.
#' @export delete.hypothesis
#' 
delete.hypothesis <- function(x, event = NA, cause = NA, effect = NA) {

    is.compliant(x)
    if (has.model(x)) {
        stop("There's a reconstructed model, hypotheses cannot be deleted now.\nUse delete.model()")
    }

    hypo_map = as.hypotheses(x)
    to_remove = c()

    if (!any(is.na(event))) {
        if (length(event) == 1 && event %in% as.events(x)[ , 'event']) {
            cause_del = which(hypo_map[ , 'cause event'] == event)
            effect_del = which(hypo_map[ , 'effect event'] == event)
            to_remove = unique(c(to_remove, cause_del, effect_del))
        } else {
            stop('Select only one event present in as.events')
        }
    }

    if (! any(is.na(cause))) {
        if (all(cause %in% as.events(x)[ , 'event'])) {
            cause_del = which(hypo_map[ , 'cause event'] == cause)
            to_remove = unique(c(to_remove, cause_del))
        } else {
            stop('Wrong cause, select only events present in as.events')
        }
    }

    if (! any(is.na(effect))) {
        if (all(effect %in% as.events(x)[ , 'event'])) {
            effect_del = which(hypo_map[ , 'effect event'] == effect)
            to_remove = unique(c(to_remove, effect_del))
        } else {
            stop('Wrong effect, select only events present in as.events')
        }
    }

    if (length(to_remove) == 0) {
        warning("Nothing to remove")
        return(x)
    }

    x$hypotheses$num.hypotheses = x$hypotheses$num.hypotheses - 1
    x$hypotheses$hlist = x$hypotheses$hlist[-to_remove, ,drop=FALSE]

    is.compliant(x)

    return(x)
}


#' Delete a pattern and every associated hypotheses from the dataset
#' @title delete.pattern
#'
#' @examples
#' data(test_dataset)
#' delete.pattern(test_dataset, pattern='XOR_EZH2')
#'
#' @param x A TRONCO compliant dataset.
#' @param pattern A pattern name
#' @return A TRONCO complian dataset.
#' @export delete.pattern
#' 
delete.pattern <- function(x, pattern) {
    is.compliant(x)
    if (has.model(x)) {
        stop("There's a reconstructed model, a pattern cannot be deleted now. \nUse delete.model()")
    }

    if (! pattern %in% as.patterns(x)) {
        stop(paste(pattern, " not in as.patterns()"))
    }

    x = delete.hypothesis(x, pattern)
    x$annotations = x$annotations[-which(rownames(x$annotations) == pattern), , drop = FALSE]
    x$genotypes = x$genotypes[, -which(colnames(x$genotypes) == pattern), drop = FALSE]
    x$hypotheses$patterns[pattern] = NULL
    x$hypotheses$pvalues = NULL

    for (atom in names(x$hypotheses$atoms)) {
        if (x$hypotheses$atoms[atom] == pattern) {
            x$hypotheses$atoms[atom] = NULL
        }
    }

    x$hypotheses$hstructure[pattern] = NULL

    if (! 'Pattern' %in% unique(x$annotations[,'type'])) {
        x$types = x$types[-which(rownames(x$types) == 'Pattern'),,drop=FALSE]

    }

    is.compliant(x)

    return(x)
}


#' Delete a reconstructed model from the dataset
#' @title delete.model
#'
#' @examples
#' data(test_model)
#' model = delete.model(test_model)
#' has.model(model)
#'
#' @param x A TRONCO compliant dataset.
#' @return A TRONCO complian dataset.
#' @export delete.model
#' 
delete.model <- function(x) {
    is.compliant(x)
    if (! has.model(x)) {
        stop("No model to delete in dataset")
    }
    x$model = NULL
    x$bootstrap = NULL
    x$confidence = NULL
    x$parameters = NULL
    x$adj.matrix.prima.facie = NULL
    x$execution.time = NULL
    x$kfold = NULL

    is.compliant(x)

    return(x)
}


#' Delete samples from selected dataset
#' @title delete.samples
#'
#' @examples
#' data(test_dataset)
#' dataset = delete.samples(test_dataset, c('patient 1', 'patient 4'))
#'
#' @param x A TRONCO compliant dataset.
#' @param samples An array of samples name
#' @return A TRONCO complian dataset.
#' @export delete.samples
#' 
delete.samples <- function(x, samples) {
    is.compliant(x, 'delete.samples input')

    if (has.model(x)) {
        stop("There's a reconstructed model, a sample cannot be deleted now. \nUse delete.model()")
    }

    samples = unique(samples)
    if (!all(samples %in% as.samples(x))) {
        stop("Not all sample in as.samples(x)")
    }

    genotypes = x$genotypes[!rownames(x$genotypes) %in% samples, , drop = FALSE]

    if (has.stages(x)) {
        stages = x$stages[!rownames(x$stages) %in% samples, , drop = FALSE]
        x$stages = stages
    }

    x$genotypes = genotypes

    is.compliant(x, 'delete.samples output')

    return(x)
}


#' Intersect samples and events of two dataset
#' @title intersect.datasets
#'
#' @examples
#' data(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @param y A TRONCO compliant dataset.
#' @param intersect.genomes If False -> just samples
#' @return A TRONCO complian dataset.
#' @export intersect.datasets
#' 
intersect.datasets <- function(x,y, intersect.genomes = TRUE) {
    is.compliant(x)
    is.compliant(y)

    ## Common samples and genes (according to intersect.genomes)
    samples = intersect(as.samples(x), as.samples(y))
    genes =
        ifelse(intersect.genomes,
               intersect(as.genes(x), as.genes(y)), # intersect.genomes -> INTERSECTION
               unique(c(as.genes(x), as.genes(y)))) # !intersect.genomes -> UNION

    report = data.frame(row.names = c('Samples', 'Genes'))
    report$x = c(nsamples(x), ngenes(x))
    report$y = c(nsamples(y), ngenes(y))

    ## Restrict genes - only if intersect.genomes = T
    
    if (intersect.genomes) {
        x = events.selection(x, filter.in.names=genes)
        y = events.selection(y, filter.in.names=genes)
    }

    ## TODO: check they have no events in common!
    ## if (as.events(x) )

    ## Restrict stamples.
    
    x = samples.selection(x, samples)
    y = samples.selection(y, samples)

    ## Result.
    
    z = ebind(x,y)

    cat('*** Intersect dataset [ intersect.genomes =', intersect.genomes, ']\n')
    report$result = c(nsamples(z), ngenes(z))
    print(report)

    return(z)
}


#' Binds events from one or more datasets, which must be defined over the same set of samples.
#' @title ebind
#'
#' @param ... the input datasets
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO complian dataset.
#' @export ebind
#' 
ebind <- function(...,
                  silent = FALSE) {
    
    ## Merge two  datasets at a time.
    events.pairwise.bind <- function(x, y) {
        is.compliant(x, 'ebind: input x')
        is.compliant(y, 'ebind: input y')

        samples.intersect = intersect(as.samples(x), as.samples(y))
        if (!(setequal(samples.intersect, as.samples(x))
              && setequal(samples.intersect, as.samples(y))))
            stop('Datasets have different samples, won\'t bind!')

        z = list()

        y$genotypes = y$genotypes[rownames(x$genotypes), , drop = FALSE]
        y$stages = y$stages[rownames(y$genotypes), , drop = FALSE]

        ## Copy genotype matrix, and sets its rownames (samples).
        
        z$genotypes = cbind(x$genotypes, y$genotypes)
        colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')

        ## Copy annotations for gene symbols etc.
        
        z$annotations = rbind(x$annotations, y$annotations)
        rownames(z$annotations) = colnames(z$genotypes)

        ## Copy types.
        
        z$types = unique(rbind(x$types, y$types))

        ## Copy stages, if present.
        
        if (has.stages(x) && has.stages(y)) {
            stages.x = as.stages(x)
            stages.y = as.stages(y)

            stages.x = stages.x[!is.na(stages.x)]
            stages.y = stages.y[!is.na(stages.y)]

            if (any(stages.x != stages.y)) {
                stop('Patients have different stages, won\'t merge!')
            }
        }

        if (has.stages(x)) {
            z = annotate.stages(z, as.stages(x))
        }

        is.compliant(z, 'ebind: output')

        return(z)
    }

    input = list(...)

    if (!silent) {
        cat('*** Binding events for', length(input), 'datasets.\n')
    }
    return(Reduce(events.pairwise.bind, input))
}


#' Binds samples from one or more datasets, which must be defined over the same set of events
#' @title sbind
#'
#' @param ... the input datasets
#' @return A TRONCO complian dataset.
#' @export sbind
#' 
sbind <- function(...) {

    ## Merge two datasets at a time.
    
    samples.pairwise.bind <- function(x, y) {
        is.compliant(x, 'sbind: input x')
        is.compliant(y, 'sbind: input y')

        if (!all(as.events(x) == as.events(y))) {
            stop('Datasets have different events, can not bind!')
        }
        
        z = list()

        ## Copy genotypes and annotations.
        
        z$genotypes = rbind(x$genotypes, y$genotypes)
        z$annotations = x$annotations
        z$types = x$types

        ## Copy stages, if present.
        
        if (has.stages(x) || has.stages(y)) {
            if (has.stages(x)) xstages = as.stages(x)
            else xstages = matrix(rep('NA', nsamples(x)), nrow = nsamples(x))

            if (has.stages(y)) ystages = as.stages(y)
            else ystages = matrix(rep('NA', nsamples(y)), nrow = nsamples(y))

            z$stages = (rbind(x$stages, y$stages))

            colnames(z$stages) = 'stage'
        }
        is.compliant(z, 'sbind: output')

        return(z)
    }

    input = list(...)
    return(Reduce(samples.pairwise.bind, input))
}


#' For an input dataset merge all the events of two or more distincit types
#' (e.g., say that missense and indel mutations are events
#' of a unique "mutation" type)
#' @title join.types
#'
#' @examples
#' data(test_dataset_no_hypos)
#' join.types(test_dataset_no_hypos, 'ins_del', 'missense_point_mutations')
#' join.types(test_dataset_no_hypos, 'ins_del',
#'      'missense_point_mutations', new.type='mut', new.color='green')
#'
#' @param x A TRONCO compliant dataset.
#' @param ... type to merge
#' @param new.type label for the new type to create
#' @param new.color color for the new type to create
#' @param silent A parameter to disable/enable verbose messages.
#' @return A TRONCO compliant dataset.
#' @export join.types
#' 
join.types <- function(x,
                       ...,
                       new.type = "new.type",
                       new.color = "khaki",
                       silent = FALSE) {

    ## Check if x is compliant.
    
    is.compliant(x)

    if(has.model(x)) {
        stop("There's a reconstructed model, types cannot be merged now. \nUse delete.model()")
    }

    if(length(as.patterns(x)) > 0) {
        stop('Patterns found. Delete patterns first.\n')
    }

    input = list(...)

    ## TODO Change this in a better way (deafult ellipsis?)
    if (length(input) == 1 && is.null(input[[1]])) {
        input = as.list(as.types(x))
    }

    if (!silent) {
        cat("*** Aggregating events of type(s) {",
            paste(unlist(input), collapse = ", ", sep = ""),
            "}\nin a unique event with label \"",
            new.type,
            "\".\n")
    }

    if (length(input) <= 1) {
        if (!silent) {
            cat("One input type provided, using renaming functions.\n")
        }

        x = rename.type(x, input[[1]], new.type)
        x = change.color(x, new.type, new.color)
        return(x)
    }


    types.check =
        lapply(input,
               function(type) {
                   type %in% as.types(x)
               })

    if (!all(unlist(types.check))) {
        t = (input[!unlist(types.check)])

        stop("No events of type '", t, "' in input dataset, will not merge.")
    }

    if (any(duplicated(input))) {
        stop("No duplicated types are allowed, will not merge.")
    }

    if (new.type %in% as.types(x)) {
        stop(paste0(new.type, "is already used in input dataset, will not merge"))
    }

    input = unlist(input)

    genes = as.genes(x, types = input)
    if (!silent) {
        cat("Dropping event types",
            paste(input, collapse = ", ", sep = ""),
            "for",
            length(genes),
            "genes.\n")
    }
    geno.matrix = matrix(, nrow = nsamples(x), ncol = length(genes))

    #if (!silent) {
    #    pb = txtProgressBar(1, length(genes), style = 3)
    #    flush.console()
    #}
    
    for (i in 1:length(genes)) {
        #if (!silent) {
        #    setTxtProgressBar(pb, i)
        #}
        if (!silent) {
            cat('.')
        }

        geno = as.matrix(rowSums(as.gene(x, genes[i], types = input)))
        geno[geno > 1] = 1

        geno.matrix[, i] = geno
    }
    if (!silent) {
        cat('\n')
    }

    rownames(geno.matrix) = as.samples(x)
    colnames(geno.matrix) = genes
    #if (!silent) {
    #    close(pb)
    #}

    z = import.genotypes(geno.matrix, event.type = new.type, color = new.color)
    if (has.stages(x)) {
        z = annotate.stages(z, as.stages(x))
    }

    y = x
    for (i in input) {
        y = delete.type(y, i)
    }

    w = ebind(y, z, silent = silent)
    is.compliant(w)

    return(w)
}


#' Deletes all events which have frequency 0 in the dataset.
#' @title trim
#'
#' @examples
#' data(test_dataset)
#' test_dataset = trim(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A TRONCO compliant dataset.
#' @export trim
#' 
trim <- function(x) {
    is.compliant(x, 'trim: input')

    x = enforce.numeric(x)

    del = names(which(colSums(x$genotypes) == 0))

    x$genotypes = x$genotypes[ , !colnames(x$genotypes) %in% del, drop = FALSE]
    x$annotations = x$annotations[!rownames(x$annotations) %in% del, , drop = FALSE]

    x$types = matrix(x$types[unique(x$annotations[ , 'type', drop = FALSE]), , drop = FALSE ], ncol = 1)
    rownames(x$types) = unique(x$annotations[,'type', drop = FALSE])
    colnames(x$types) = 'color'

    is.compliant(x, 'trim: output')
    return(x)
}


# TODO: check
#' Split cohort (samples) into groups, return either all groups or a specific group.
#' @title ssplit
#'
#' @param x A TRONCO compliant dataset.
#' @param clusters A list of clusters. Rownames must match samples list of x
#' @param idx ID of a specific group present in stages. If NA all groups will be extracted
#' @return A TRONCO compliant dataset.
#' @export ssplit
#' 
ssplit <- function(x, clusters, idx = NA) {
    is.compliant(x)

    data = x$genotypes

    cat('*** Splitting cohort into groups.\n')

    ## TODO: check that clusters has at least 2 columns

    ## Check that map has correct size.

    if (nsamples(x) != nrow(clusters))
        stop(paste("Error: cannot split, number of samples (",
                   nsamples(x) ,
                   ") and groups (",
                   nrow(clusters),
                   ") do not match.",
                   sep = ''));

    ## Check it is an actual map.
    
    if (!all(rownames(clusters) %in% rownames(data)))
        stop(paste('Error: samples',
                   paste(rownames(clusters)[!rownames(clusters) %in% as.samples(x)],
                         collapse=', ',
                         sep = ''),
                   'are not assigned to a group',
                   sep = ''));

    ## Groups info.
    
    cluster.labels = unique(clusters)
    num.clusters = nrow(cluster.labels)

    ## Extract a specific group.
    
    if (!is.na(idx)) {
        y = list()
        samples.in.cluster = rownames(clusters)[clusters == idx]

        cat(paste('Group \"', idx, '\" has ', length(samples.in.cluster),
                  ' samples, returning this group.\n', sep=''))

        y$genotypes = data[samples.in.cluster, ];
        y$annotations = x$annotations
        y$types = x$types

        if (!is.null(x$stages)) {
            y$stages = as.matrix(x$stages[samples.in.cluster, ])
            rownames(y$stages) = samples.in.cluster
        }

        is.compliant(y, 'ssplit.split with index')
        return(y)
    }

    ## Extract all groups.
    
    partitions = list()
    for (i in 1:num.clusters) {
        y = list()
        samples.in.cluster = rownames(clusters)[clusters == cluster.labels[i, 1]]

        cat(paste('Group \"', cluster.labels[i,1], '\" has ', length(samples.in.cluster),
                  ' samples.\n', sep=''))

        y$genotypes = data[samples.in.cluster, ];
        y$annotations = x$annotations
        y$types = x$types

        if (!is.null(x$stages)) {
            y$stages = as.matrix(x$stages[samples.in.cluster, ])
            rownames(y$stages) = samples.in.cluster
        }

        is.compliant(y, 'subtypes.split partitionig')
        partitions = append(partitions, list(y))
        names(partitions)[i] = cluster.labels[i,1]
    }

    return(partitions)
}


# Split events into groups according to their types.
#
# x: cohort
# @export

tsplit <- function(x) {
    ## Parametro per estrarre un tipo solo di evento?
}


#' Merge a list of events in an unique event
#' @title join.events
#'
#' @examples
#' data(muts)
#' dataset = join.events(muts, 'G1', 'G2', new.event='test', new.type='banana', event.color='yellow')
#'
#' @param x A TRONCO compliant dataset.
#' @param ... A list of events to merge
#' @param new.event The name of the resultant event
#' @param new.type The type of the new event
#' @param event.color The color of the new event
#' @return A TRONCO compliant dataset.
#' @export join.events
#' 
join.events <- function(x, ..., new.event, new.type, event.color) {

    events = list(...)

    if (length(events) < 2) {
        stop('Badly formed events.')
    }

    for (pattern in as.patterns(x)) {
        if (any(events %in% rownames(as.events.in.patterns(x, patterns=pattern)))) {
            stop('Found event in pattern \"', pattern, '\". Delete that pattern first.\n')
        }
    }

    x = enforce.numeric(x)

    y = x
    as_ev = as.events(x)

    for (event in events) {
        y = delete.event(y, gene = as_ev[event, 'event'], type = as_ev[event, 'type'])
    }

    genos = x$genotypes[, unlist(events)]
    genos = matrix(rowSums(genos), nrow = nsamples(x))
    genos[genos > 1] = 1
    colnames(genos) = new.event
    rownames(genos) = as.samples(x)

    genos = import.genotypes(genos, event.type = new.type, color = event.color)
    y = ebind(y, genos)

    return(y)
}


# rename field of reconstruction
rename.reconstruction.fields <- function(reconstruction, genotypes){
    rownames(reconstruction$confidence) =
        c("temporal priority",
          "probability raising",
          "hypergeometric test")
    colnames(reconstruction$confidence) = "confidence"
    rownames(reconstruction$confidence[[1,1]]) = colnames(genotypes)
    colnames(reconstruction$confidence[[1,1]]) = colnames(genotypes)
    rownames(reconstruction$confidence[[2,1]]) = colnames(genotypes)
    colnames(reconstruction$confidence[[2,1]]) = colnames(genotypes)
    rownames(reconstruction$confidence[[3,1]]) = colnames(genotypes)
    colnames(reconstruction$confidence[[3,1]]) = colnames(genotypes)

    for (i in 1:length(reconstruction$model)) {

        ## Set rownames and colnames to the probabilities.
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = colnames(genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$marginal.probs) = "marginal probability"
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$joint.probs) = colnames(genotypes)
        rownames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = colnames(genotypes)
        colnames(reconstruction$model[[i]]$probabilities$probabilities.observed$conditional.probs) = "conditional probability"

        ## Set rownames and colnames to the parents positions.
        rownames(reconstruction$model[[i]]$parents.pos) = colnames(genotypes)
        colnames(reconstruction$model[[i]]$parents.pos) = "parents"

        ## Set rownames and colnames to the adjacency matrices.
        rownames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(genotypes)
        colnames(reconstruction$model[[i]]$adj.matrix$adj.matrix.fit) = colnames(genotypes)

    }
    return(reconstruction)
}


#### end of file -- editing.functions.R
