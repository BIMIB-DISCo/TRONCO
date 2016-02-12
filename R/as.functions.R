#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Return all genotypes for input 'x', which should be a TRONCO compliant dataset
#' see \code{is.compliant}. 
#' Function \code{keysToNames} can be used to translate colnames to events. 
#' @title as.genotypes
#'
#' @examples
#' data(test_dataset)
#' as.genotypes(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A TRONCO genotypes matrix.
#' @export as.genotypes
#' 
as.genotypes <- function(x) {
    return(x$genotypes)
}


#' Return all sample IDs for input 'x', which should be a TRONCO compliant dataset - see \code{is.compliant}. 
#' @title as.samples
#'
#' @examples
#' data(test_dataset)
#' as.samples(test_dataset)
#' 
#' @param x A TRONCO compliant dataset.
#' @return A vector of sample IDs
#' @export as.samples
#' 
as.samples <- function(x) {
    return(rownames(x$genotypes))
}


#' Return all gene symbols for which a certain type of event exists in 'x', which should be a
#' TRONCO compliant dataset - see \code{is.compliant}. 
#' @title as.genes
#'
#' @examples
#' data(test_dataset)
#' as.genes(test_dataset)
#' 
#' @param x A TRONCO compliant dataset.
#' @param types The types of events to consider, if NA all available types are used.
#' @return A vector of gene symbols for which a certain type of event exists
#' @export as.genes
#' 
as.genes <- function(x, types = NA) {
    if ('Pattern' %in% types) {
        stop('"Pattern" is not a valid gene type, it is a reseverd keyword in TRONCO.')
    }
    
    if (npatterns(x) > 0) {
        ev = as.events(x, types = types)
        ev = ev[which(ev[ , 'type'] != 'Pattern'), 'event']
        return(unique(ev))
    } else
        return(unique(as.events(x, types = types)[ , 'event']))
}


#' Return all events involving certain genes and of a certain type in 'x', which should be a
#' TRONCO compliant dataset - see \code{is.compliant}. 
#' @title as.events
#'
#' @examples
#' data(test_dataset)
#' as.events(test_dataset)
#' as.events(test_dataset, types='ins_del')
#' as.events(test_dataset, genes = 'TET2')
#' as.events(test_dataset, types='Missing')
#'
#' @param x A TRONCO compliant dataset.
#' @param types The types of events to consider, if NA all available types are used.
#' @param genes The genes to consider, if NA all available genes are used.
#' @param keysToNames If TRUE return a list of mnemonic name composed by type + gene
#' @return A matrix with 2 columns (event type, gene name) for the events found.
#' @export as.events
#' 
as.events <- function(x, genes = NA, types = NA, keysToNames = FALSE) {
    ann = x$annotations[ , c('type', 'event'), drop = FALSE]

    if (!any(is.na(genes)))
        ann =
            ann[which(ann[, 'event', drop=FALSE] %in% genes), , drop = FALSE]
    
    if (!any(is.na(types)))
        ann =
            ann[which(ann[, 'type', drop = FALSE] %in% types), , drop = FALSE] 

    if (keysToNames) {
        ann = keysToNames(x, ann)
        names(colnames(ann)) = NULL
        names(rownames(ann)) = NULL
        return(rownames(ann))
    }  

    return(ann)
}


#' Return the association sample -> stage, if any. Input 'x' should be a
#' TRONCO compliant dataset - see \code{is.compliant}.
#' @title as.stages
#'
#' @examples
#' data(test_dataset)
#' data(stage)
#' test_dataset = annotate.stages(test_dataset, stage)
#' as.stages(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A matrix with 1 column annotating stages and rownames as sample IDs.
#' @export as.stages
#' 
as.stages <- function(x) {
    if (has.stages(x)) {
        return(x$stages) 
    }
    return(NA)
}


#' Return the types of events for a set of genes which are in 'x', which should be a
#' TRONCO compliant dataset - see \code{is.compliant}.
#' @title as.types
#'
#' @examples
#' data(test_dataset)
#' as.types(test_dataset)
#' as.types(test_dataset, genes='TET2')
#' 
#' @param x A TRONCO compliant dataset.
#' @param genes A list of genes to consider, if NA all genes are used.
#' @return A matrix with 1 column annotating stages and rownames as sample IDs.
#' @export as.types
#' 
as.types <- function(x, genes = NA) {    
    return(unlist(unique(as.events(x, genes = genes)[, 'type'])))
}


#' Return the colors associated to each type of event in 'x', which should be a
#' TRONCO compliant dataset - see \code{is.compliant}.
#' @title as.colors
#' 
#' @examples
#' data(test_dataset)
#' as.colors(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return A named vector of colors.
#' @export as.colors
#' 
as.colors <- function(x) {
    return(x$types[, 'color'])
}


#' Return the genotypes for a certain set of genes and type of events. Input 'x' should be a
#' TRONCO compliant dataset - see \code{is.compliant}. In this case column names are substituted
#' with events' types.
#' @title as.gene
#'
#' @examples
#' data(test_dataset)
#' as.gene(test_dataset, genes = c('EZH2', 'ASXL1'))
#' 
#' @param x A TRONCO compliant dataset.
#' @param types The types of events to consider, if NA all available types are used.
#' @param genes The genes to consider, if NA all available genes are used.
#' @return A matrix, subset of \code{as.genotypes(x)} with colnames substituted  with events' types.
#' @export as.gene
#' 
as.gene <- function(x, genes, types = NA) {
    keys = as.events(x, genes = genes, types = types)

    data =
        data.frame(x$genotypes[, rownames(keys)],
                   row.names = as.samples(x))
    
    colnames(data) =
        apply(keys, 1, FUN = paste, collapse = ' ')

    return(data)
}


#' Return a dataset where all events for a gene are merged in a unique event, i.e., 
#' a total of gene-level alterations diregarding the event type. Input 'x' is checked
#' to be a TRONCO compliant dataset - see \code{is.compliant}. 
#' @title as.alterations
#'
#' @examples
#' data(muts)
#' as.alterations(muts)
#' 
#' @param x A TRONCO compliant dataset.
#' @param new.type The types label of the new event type, 'Alteration' by default.
#' @param new.color The color of the event \code{new.type}, default 'khaki'.
#' @return A TRONCO compliant dataset with alteration profiles.
#' @export as.alterations
#' 
as.alterations <- function(x, new.type = 'Alteration', new.color = 'khaki') {
    is.compliant(x)
    if(has.model(x)) {
        stop("There's a reconstructed model, types cannot be merged now. \nUse delete.model()")
    }

    if(length(as.patterns(x)) > 0) {
        stop('Patterns found. Delete patterns first.\n')
    }

    merge.types(x, NULL, new.type = new.type, new.color = new.color)
}

#' Return the patterns in the dataset which constitute CAPRI's hypotheses.
#' @title as.patterns
#'
#' @examples
#' data(test_dataset)
#' as.patterns(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @return The patterns in the dataset which constitute CAPRI's hypotheses.
#' @export as.patterns
#' 
as.patterns <- function(x) {
    if (length(x$hypotheses) == 0 || is.na(x$hypotheses)) {
        return(NULL)
    }

    is.compliant(x)
    if ('hstructure' %in% names(x$hypotheses)) {
        return(ls(x$hypotheses$hstructure))
    }
    return(NULL)
}


#' Return the hypotheses in the dataset which constitute CAPRI's hypotheses.
#' @title as.hypotheses
#'
#' @examples
#' data(test_dataset)
#' as.hypotheses(test_dataset)
#'
#' @param x A TRONCO compliant dataset.
#' @param cause A list of genes to use as causes
#' @param effect A list of genes to use as effects
#' @return The hypotheses in the dataset which constitute CAPRI's hypotheses.
#' @export as.hypotheses
#' 
as.hypotheses <- function(x, cause = NA, effect = NA) {
    if (nhypotheses(x) < 1)
        return(NULL)

    hlist = x$hypotheses$hlist

    list_c = x$annotations[hlist[,'cause'], c('type', 'event'), drop = F]
    colnames(list_c) = c('cause type', 'cause event')
    rownames(list_c) = NULL
    list_e = x$annotations[hlist[,'effect'], c('type', 'event'), drop = F]
    colnames(list_e) = c('effect type', 'effect event')
    rownames(list_e) = NULL

    filtered_list = cbind(list_c, list_e)

    if (!is.na(cause)) {
        if (all(cause %in% as.events(x)[, 'event'])) {
            filtered_list = filtered_list[which(filtered_list[, 'cause event'] == cause), ]
        } else {
            stop('some cause not in as.events\n')
        }
    }

    if (!is.na(effect)) {
        if (all(effect %in% as.events(x)[, 'event'])) {
            filtered_list = filtered_list[which(filtered_list[, 'effect event'] == effect), ]
        } else {
            stop('some effect not in as.events\n')
        }
    }

    return(filtered_list)
}


#' Return the list of events present in selected patterns
#'
#' @examples
#' data(test_dataset)
#' as.events.in.patterns(test_dataset)
#' as.events.in.patterns(test_dataset, patterns='XOR_EZH2')
#'
#' @title as.events.in.patterns
#' @param x A TRONCO compliant dataset.
#' @param patterns A list of patterns for which the list will be returned
#' @return A list of events present in patterns which consitute CAPRI's hypotheses
#' @export as.events.in.patterns
#' 
as.events.in.patterns <- function(x, patterns = NULL) {
    is.compliant(x)
    ann = x$annotations[, c('type', 'event'), drop = FALSE]
    if (is.null(patterns)) {
        patterns = as.patterns(x)
    }

    genes_list = NULL
    for (h in patterns) {
        if (!h %in% as.patterns(x)) {
            stop('Hypothesis ', h, ' not in as.patterns(x)')
        }

        g =
            lapply(colnames(x$hypotheses$hstructure[[h]]),
                   function(x) { if (!is.logic.node(x)) return(x) })
        genes_list = append(genes_list, g)
    }  
    genes_list = unique(unlist(genes_list))

    if (!(is.null(genes_list)))
        ann = ann[ which(rownames(ann) %in% genes_list), , drop = FALSE] 

    return(ann)
}


#' Return the list of genes present in selected patterns
#'
#' @examples
#' data(test_dataset)
#' as.genes.in.patterns(test_dataset)
#' as.genes.in.patterns(test_dataset, patterns='XOR_EZH2')
#'
#' @title as.genes.in.patterns
#' @param x A TRONCO compliant dataset.
#' @param patterns A list of patterns for which the list will be returned
#' @return A list of genes present in patterns which consitute CAPRI's hypotheses
#' @export as.genes.in.patterns
#' 
as.genes.in.patterns <- function(x, patterns=NULL) {
    events = as.events.in.patterns(x, patterns)
    genes = unique(events[, 'event'])
    return(genes)
}


#' Return the list of types present in selected patterns
#'
#' @examples
#' data(test_dataset)
#' as.types.in.patterns(test_dataset)
#' as.types.in.patterns(test_dataset, patterns='XOR_EZH2')
#'
#' @title as.types.in.patterns
#' @param x A TRONCO compliant dataset.
#' @param patterns A list of patterns for which the list will be returned
#' @return A list of types present in patterns which consitute CAPRI's hypotheses
#' @export as.types.in.patterns
#' 
as.types.in.patterns <- function(x, patterns = NULL) {
    events = as.events.in.patterns(x, patterns)
    types = unique(events[, 'type'])
    return(types)
}


#' Return a list of events which are observed in the input samples list
#'
#' @examples
#' data(test_dataset)
#' as.events.in.sample(test_dataset, c('patient 1', 'patient 7'))
#'
#' @title as.events.in.sample
#' @param x A TRONCO compliant dataset
#' @param sample Vector of sample names
#' @return A list of events which are observed in the input samples list
#' @export as.events.in.sample
#' 
as.events.in.sample <- function(x, sample) {
    
    aux <- function(s) {
        sub.geno = as.genotypes(x)[s, , drop = FALSE]  
        sub.geno = sub.geno[, sub.geno == 1, drop = FALSE]
        return(as.events(x)[colnames(sub.geno), , drop = FALSE])
    }

    return(sapply(sample, FUN = aux))
}


#' Return confidence information for a TRONCO model. Available information are: temporal priority (tp), 
#' probability raising (pr), hypergeometric test (hg), parametric (pb), non parametric (npb) or 
#' statistical (sb) bootstrap, entropy loss (eloss), prediction error (prederr).
#' Confidence is available only once a model has been reconstructed with any of the algorithms implemented
#' in TRONCO. If more than one model has been reconstructed - for instance via multiple regularizations - 
#' confidence information is appropriately nested. The requested confidence is specified via 
#' vector parameter \code{conf}.
#'
#' @examples
#' data(test_model)
#' as.confidence(test_model, conf='tp')
#' as.confidence(test_model, conf=c('tp', 'hg'))
#'
#' @title as.confidence
#' @param x A TRONCO model.
#' @param conf A vector with any of 'tp', 'pr', 'hg', 'npb', 'pb', 'sb', 'eloss', 'prederr' or 'posterr'. 
#' @return A list of matrices with the event-to-event confidence. 
#' @export as.confidence
#' 
as.confidence <- function(x, conf) {
    is.compliant(x)
    is.model(x)
    if (!is.vector(conf)) stop('"conf" should be a vector.')

    keys = c('hg', 'tp', 'pr', 'npb', 'pb', 'sb', 'eloss', 'prederr', 'posterr')

    if (!all(conf %in% keys)) 
        stop('Confidence keyword unrecognized, \'conf\' should be any of:\n
            INPUT DATASET\n
            \t \"hg\" - hypergeometric test (randomness of observations)\n
            SELECTIVE ADVANTAGE SCORES\n
            \t \"tp\" - temporal priority (temporal ordering of events)\n
            \t \"pr\" - probability raising (selectivity among events)\n
            MODEL CONFIDENCE - requires post-reconstruction bootstrap or cross-validation\n
            \t \"npb\"     - non-parametric bootstrap,\n
            \t \"pb\"      - parametric bootstrap\n
            \t \"sb\"      - statistical bootstrap\n
            \t \"eloss\"   - entropy loss\n
            \t \"prederr\" - prediction error\n
            \t \"postderr\"- posterior classification error'
             )

    if (is.null(x$confidence)
        || is.na(x$confidence)
        || is.null(x$model)
        || is.na(x$model))
        stop('Input \'x\' does not contain a TRONCO model. No confidence to show.\n')

    models = names(x$model)

#     if (is.null(x$bootstrap))
#       stop('No bootstrap executed in this TRONCO object.')
#         
    has.npb.bootstrap = is.null(x$bootstrap[[models[1]]]$npb)
    has.pb.bootstrap = is.null(x$bootstrap[[models[1]]]$pb)
    has.sb.bootstrap = is.null(x$bootstrap[[models[1]]]$sb)

    if ('npb' %in% conf && has.npb.bootstrap) {
        stop('Non-parametric bootstrap was not performed. Remove keyword\n')
    }

    if ('pb' %in% conf && has.pb.bootstrap) {
        stop('Parametric bootstrap was not performed. Remove keyword\n')
    }

    if ('sb' %in% conf && has.sb.bootstrap) {
        stop('Statistical bootstrap was not performed. Remove keyword\n')
    }
    
    result = NULL

    if ('hg' %in% conf) {
        result$hg = x$confidence['hypergeometric test', ][[1]]
    }
    
    if ('tp' %in% conf) {
        result$tp = x$confidence['temporal priority', ][[1]]
    }
    
    if ('pr' %in% conf) {
        result$pr = x$confidence['probability raising', ][[1]]
    }
    
    if ('npb' %in% conf) {
        for (i in 1:length(models)) {
            result$npb[models[i]] =
                list(x$bootstrap[[models[i]]]$npb$bootstrap.edge.confidence)  
        }
    }
    
    if ('pb' %in% conf) {
        for (i in 1:length(models)) {
            result$pb[models[i]] =
                list(x$bootstrap[[models[i]]]$pb$bootstrap.edge.confidence)  
        }
    }

    if ('sb' %in% conf) { 
        for (i in 1:length(models)) {
            result$sb[models[i]] =
                list(x$bootstrap[[models[i]]]$sb$bootstrap.edge.confidence)  
        }
    }
    
    if ('eloss' %in% conf) { 
        result$eloss = as.kfold.eloss(x, models)
    }

    if ('prederr' %in% conf) {
        result$prederr = as.kfold.prederr(x, models = models)
    }

    if ('posterr' %in% conf) {
        result$posterr = as.kfold.posterr(x, models = models)
    }
    
    return(result)  
}


#' Extract the models from a reconstructed object.
#'
#' @examples
#' data(test_model)
#' as.models(test_model)
#'
#' @title as.models
#' @param x A TRONCO model.
#' @param models The name of the models to extract, e.g. 'bic', 'aic', 'caprese', all by default. 
#' @return The models in a reconstructed object. 
#' @export as.models
#' 
as.models <- function(x, models=names(x$model)) {
    is.compliant(x)
    is.model(x)
    if (!is.vector(models)) {
        stop('"models" should be a vector.')
    }

    for (model in models) {
        if ( !model %in% names(x$model)) {
            stop(paste('model:', model, 'not present'))
        }
    }

    return(x$model[models])
}


#' Return the description annotating the dataset, if any. Input 'x' should be
#' a TRONCO compliant dataset - see \code{is.compliant}. 
#'
#' @examples
#' data(test_dataset)
#' as.description(test_dataset)
#'
#' @title as.description
#' @param x A TRONCO compliant dataset.
#' @return The description annotating the dataset, if any.
#' @export as.description
#' 
as.description <- function(x) {
    if (!is.null(x$name))
        return(x$name)
    return("")
}


#' Given a cohort and a pathway, return the cohort with events restricted to genes 
#' involved in the pathway. This might contain a new 'pathway' genotype with an alteration mark if
#' any of the involved genes are altered. 
#'
#' @examples
#' data(test_dataset)
#' p = as.pathway(test_dataset, c('ASXL1', 'TET2'), 'test_pathway')
#'
#' @title as.pathway
#' @param x A TRONCO compliant dataset.
#' @param pathway.genes Gene (symbols) involved in the pathway.
#' @param pathway.name Pathway name for visualization.
#' @param pathway.color Pathway color for visualization.
#' @param aggregate.pathway If TRUE drop the events for the genes in the pathway.
#' @return Extract the subset of events for genes which are part of a pathway.
#' @export as.pathway
#' 
as.pathway <- function(x, 
                       pathway.genes, 
                       pathway.name, 
                       pathway.color='yellow', 
                       aggregate.pathway = TRUE) {
    
    is.compliant(x, 'as.pathway: input')
    data = x$genotypes
    cat(paste('*** Extracting events for pathway: ',
              pathway.name,
              '.\n',
              sep = ''))

    ## Select only those events involving a gene in pathway.genes
    ## which is also in x.
    
    y = events.selection(x, NA, filter.in.names = pathway.genes, NA)

    ## Extend genotypes.
    
    y = enforce.numeric(y)

    pathway =
        data.frame(rowSums(as.genotypes(y)),
                   row.names = as.samples(y),
                   stringsAsFactors = FALSE)  
    pathway[pathway > 1, ] =  1
    colnames(pathway) = pathway.name 

    pathway =
        import.genotypes(pathway,
                         event.type = 'Pathway',
                         color = pathway.color)

    cat('Pathway extracted succesfully.\n')

    if (!aggregate.pathway) {
        pathway = ebind(pathway, y)
    }

    if (has.stages(y)) {
        pathway = annotate.stages(pathway, as.stages(y))
    }

    is.compliant(pathway, 'as.pathway: output')

    return(pathway)
}


#' Extract the adjacency matrix of a TRONCO model. The matrix is indexed with colnames/rownames which 
#' represent genotype keys - these can be resolved with function \code{keysToNames}. It is possible to
#' specify a subset of events to build the matrix, a subset of models if multiple reconstruction have
#' been performed. Also, either the prima facie matrix or the post-regularization matrix can be extracted.
#'
#' @examples
#' data(test_model)
#' as.adj.matrix(test_model)
#' as.adj.matrix(test_model, events=as.events(test_model)[5:15,])
#' as.adj.matrix(test_model, events=as.events(test_model)[5:15,], type='pf')
#' 
#' @title as.adj.matrix
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param type Either the prima facie ('pf') or the post-regularization ('fit') matrix, 'fit' by default.
#' @return The adjacency matrix of a TRONCO model. 
#' @export as.adj.matrix
#' 
as.adj.matrix <- function(x,
                          events = as.events(x),
                          models = names(x$model),
                          type = 'fit') {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    if (!is.vector(models)) {
        stop('"models" should be a vector.') 
    }
    
    if (!type %in% c('fit', 'pf')  ) {
        stop('"type" should be any of \'fit\' (post-regularization) or \'pf\' (prima facie).')
    }  

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        if (type == 'fit') mat = m[[i]]$adj.matrix$adj.matrix.fit
        if (type == 'pf') mat = m[[i]]$adj.matrix$adj.matrix.pf

        mat = mat[rownames(events), , drop = FALSE]
        mat = mat[, rownames(events), drop = FALSE]

        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


#' Extract the marginal probabilities from a TRONCO model. The return matrix is indexed with rownames which 
#' represent genotype keys - these can be resolved with function \code{keysToNames}. It is possible to
#' specify a subset of events to build the matrix, a subset of models if multiple reconstruction have
#' been performed. Also, either the observed or fit probabilities can be extracted.
#'
#' @examples
#' data(test_model)
#' as.marginal.probs(test_model)
#' as.marginal.probs(test_model, events=as.events(test_model)[5:15,])
#'
#' @title as.marginal.probs
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param type Either observed ('observed') or fit ('fit') probabilities, 'observed' by default.
#' @return The marginal probabilities in a TRONCO model. 
#' @export as.marginal.probs
as.marginal.probs <- function(x,
                              events = as.events(x),
                              models = names(x$model),
                              type = 'observed') {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    if (!type %in% c('observed', 'fit')) {
        stop('Marginal probabilities are available for \'observed\' (empirical) or \'fit\' (estimated).') 
    }
    
    if (any(is.null(colnames(events)))) {
        stop('Events should have rownames to access the adjacency matrix - use \'as.events\' function?')
    }

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        if (type == 'observed')
            mat =
                m[[i]]$probabilities$probabilities.observed$marginal.probs
        
        if (type == 'fit')
            mat =
                m[[i]]$probabilities$probabilities.fit$estimated.marginal.probs

        if (type == 'fit' && is.na(mat))
            stop('Marginal probabilities have not been estimated yet - see TRONCO Manual.')       

        mat = mat[rownames(events), , drop = FALSE]
        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


#' Extract the joint probabilities from a TRONCO model. The return matrix is indexed with rownames/colnames which 
#' represent genotype keys - these can be resolved with function \code{keysToNames}. It is possible to
#' specify a subset of events to build the matrix, a subset of models if multiple reconstruction have
#' been performed. Also, either the observed or fit probabilities can be extracted.
#'
#' @examples
#' data(test_model)
#' as.joint.probs(test_model)
#' as.joint.probs(test_model, events=as.events(test_model)[5:15,])
#'
#' @title as.joint.probs
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param type Either observed ('observed') or fit ('fit') probabilities, 'observed' by default.
#' @return The joint probabilities in a TRONCO model. 
#' @export as.joint.probs
#'
as.joint.probs <- function(x,
                           events = as.events(x),
                           models = names(x$model),
                           type = 'observed') {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    if (!type %in% c('observed', 'fit')) {
        stop('Joint probabilities are available for \'observed\' (empirical) or \'fit\' (estimated).') 
    }
    
    if (any(is.null(colnames(events)))) {
        stop('Events should have rownames to access the adjacency matrix - use \'as.events\' function?')
    }

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        if (type == 'observed')
            mat =
                m[[i]]$probabilities$probabilities.observed$joint.probs
        
        if (type == 'fit')
            mat =
                m[[i]]$probabilities$probabilities.fit$estimated.joint.probs

        if (type == 'fit' && is.na(mat))
            stop('Joint probabilities have not been estimated yet - see the TRONCO Manual.')        

        mat = mat[rownames(events), , drop = FALSE]
        mat = mat[, rownames(events), drop = FALSE]
        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


#' Extract the conditional probabilities from a TRONCO model. The return matrix is indexed with rownames which 
#' represent genotype keys - these can be resolved with function \code{keysToNames}. It is possible to
#' specify a subset of events to build the matrix, a subset of models if multiple reconstruction have
#' been performed. Also, either the observed or fit probabilities can be extracted.
#'
#' @title as.conditional.probs
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param type Either observed ('observed') or fit ('fit') probabilities, 'observed' by default.
#' @return The conditional probabilities in a TRONCO model. 
#' @export as.conditional.probs
#' 
as.conditional.probs <- function(x,
                                 events = as.events(x),
                                 models = names(x$model),
                                 type = 'observed') {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    if (!type %in% c('observed', 'fit'))
        stop('Conditional probabilities are available for \'observed\' (empirical) or \'fit\' (estimated).')  
    
    if (any(is.null(colnames(events))))
        stop('Events should have rownames to access the adjacency matrix - use \'as.events\' function?')

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        if (type == 'observed')
            mat =
                m[[i]]$probabilities$probabilities.observed$conditional.probs
        
        if (type == 'fit')
            mat =
                m[[i]]$probabilities$probabilities.fit$estimated.conditional.probs

        if (type == 'fit' && is.na(mat))
            stop('Conditional probabilities have not been estimated yet - see the TRONCO Manual.')        

        mat = mat[rownames(events), , drop = FALSE]
        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


# Extract the estimated rates of false positives an negatives in the data, given the model. 
# A subset of models if multiple reconstruction have been performed can be extracted.
# 
# @examples
# data(test_model)
# as.error.rates(test_model)
#
# @title as.error.rates
# @param x A TRONCO model.
# @param models A subset of reconstructed models, all by default.
# @return The estimated rates of false positives an negatives in the data, given the model.
#
as.error.rates <- function(x, models = names(x$model)) {
    is.compliant(x)
    is.model(x)

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        mat = m[[i]]$error.rates
        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


#' Returns a dataframe with all the selective advantage relations in a 
#' TRONCO model. Confidence is also shown - see \code{as.confidence}. It is possible to
#' specify a subset of events or models if multiple reconstruction have
#' been performed. 
#'
#' @examples
#' data(test_model)
#' as.selective.advantage.relations(test_model)
#' as.selective.advantage.relations(test_model, events=as.events(test_model)[5:15,])
#' as.selective.advantage.relations(test_model, events=as.events(test_model)[5:15,], type='pf')
#'
#' @title as.selective.advantage.relations
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param type Either Prima Facie ('pf') or fit ('fit') probabilities, 'fit' by default.
#' @return All the selective advantage relations in a TRONCO model 
#' @export as.selective.advantage.relations
#' 
as.selective.advantage.relations <- function(x,
                                             events = as.events(x),
                                             models = names(x$model),
                                             type = 'fit') {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    ## TEMPORARY HORRIBLE FIX
    matrix = NULL
    
    if (type == 'pf')
        matrix$pf = x$adj.matrix.prima.facie   
    else
        matrix =
            as.adj.matrix(x,
                          events = events,
                          models = models,
                          type = type)
    
    matrix = lapply(matrix, keysToNames, x = x)

    conf = as.confidence(x, conf = c('tp', 'pr', 'hg'))
    conf = lapply(conf, keysToNames, x = x)

    matrix.to.df <- function(m) {    
        entries = length(which(m == 1))
        df = NULL
        df$SELECTS = NULL
        df$SELECTED = NULL
        df$OBS.SELECTS = NULL
        df$OBS.SELECTED = NULL
        df$HG = NULL
        df$TP = NULL
        df$PR = NULL
        
        if (entries == 0) {
            return(NULL)
        }

        for (i in 1:ncol(m)) {
            for (j in 1:nrow(m)) {
                if (m[i, j] == 1) {
                    df$SELECTS = c(df$SELECTS, rownames(m)[i])
                    df$SELECTED = c(df$SELECTED, colnames(m)[j])

                    df$OBS.SELECTS =
                        c(df$OBS.SELECTS,
                          sum(as.genotypes(x)[, nameToKey(x, rownames(m)[i])]))
                    df$OBS.SELECTED =
                        c(df$OBS.SELECTED,
                          sum(as.genotypes(x)[, nameToKey(x, colnames(m)[j])]))

                    df$TP = c(df$TP, conf$tp[rownames(m)[i], colnames(m)[j]])
                    df$PR = c(df$PR, conf$pr[rownames(m)[i], colnames(m)[j]])
                    df$HG = c(df$HG, conf$hg[rownames(m)[i], colnames(m)[j]])   
                }
            }
        }

        df =
            cbind(df$SELECTS,
                  df$SELECTED,
                  df$OBS.SELECTS,
                  df$OBS.SELECTED,
                  df$TP,
                  df$PR,
                  df$HG)

        colnames(df) =
            c('SELECTS',
              'SELECTED',
              'OBS.SELECTS',
              'OBS.SELECTED',
              'TEMPORAL.PRIORITY',
              'PROBABILITY.RAISING',
              'HYPERGEOMETRIC')
        
        rownames(df) = paste(1:nrow(df))

        df = data.frame(df, stringsAsFactors = FALSE) 
        df$OBS.SELECTS = as.numeric(df$OBS.SELECTS)
        df$OBS.SELECTED = as.numeric(df$OBS.SELECTED)
        df$HYPERGEOMETRIC = as.numeric(df$HYPERGEOMETRIC)
        df$TEMPORAL.PRIORITY = as.numeric(df$TEMPORAL.PRIORITY)
        df$PROBABILITY.RAISING = as.numeric(df$PROBABILITY.RAISING)

        return(df)
    }
    
    return(lapply(matrix, matrix.to.df))
}


#' Returns a dataframe with all the bootstrap score in a 
#' TRONCO model. It is possible to specify a subset of events
#' or models if multiple reconstruction have been performed. 
#'
#' @examples
#' data(test_model)
#' as.bootstrap.scores(test_model)
#' as.bootstrap.scores(test_model, events=as.events(test_model)[5:15,])
#'
#' @title as.bootstrap.scores
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @return All the bootstrap scores in a TRONCO model 
#' @export as.bootstrap.scores
#' 
as.bootstrap.scores <- function(x,
                                events = as.events(x),
                                models = names(x$model)) {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    matrix = 
    as.adj.matrix(x,
        events = events,
        models = models,
        type = 'fit')


    matrix = lapply(matrix, keysToNames, x = x)

    if (!(is.null(x$confidence) && is.na(x$confidence))) {
        models = names(x$model)
    }

    has.npb.bootstrap = !is.null(x$bootstrap[[models[1]]]$npb)
    has.pb.bootstrap = !is.null(x$bootstrap[[models[1]]]$pb)
    has.sb.bootstrap = !is.null(x$bootstrap[[models[1]]]$sb)

    if(has.npb.bootstrap) {
        npb.boot.conf = lapply(as.confidence(x, conf = c('npb'))$npb, keysToNames, x = x)
    }
    if(has.pb.bootstrap) {
        pb.boot.conf = lapply(as.confidence(x, conf = c('pb'))$pb, keysToNames, x = x)
    }
    if(has.sb.bootstrap) {
        sb.boot.conf = lapply(as.confidence(x, conf = c('sb'))$sb, keysToNames, x = x)
    }

    matrix.to.df <- function(z) {   
        m = matrix[[z]]

        entries = length(which(m == 1))
        df = NULL
        df$SELECTS = NULL
        df$SELECTED = NULL
        df$OBS.SELECTS = NULL
        df$OBS.SELECTED = NULL
        df$BOOT.NPB = NULL
        df$BOOT.PB = NULL
        df$BOOT.SB = NULL

        if (entries == 0) {
            return(NULL)
        }

        for (i in 1:ncol(m)) {
            for (j in 1:nrow(m)) {

                if (m[i, j] == 1) {
                    df$SELECTS = c(df$SELECTS, rownames(m)[i])
                    df$SELECTED = c(df$SELECTED, colnames(m)[j])

                    df$OBS.SELECTS = 
                        c(df$OBS.SELECTS,
                          sum(as.genotypes(x)[, nameToKey(x, rownames(m)[i])])
                          )
                    
                    df$OBS.SELECTED =
                        c(df$OBS.SELECTED,
                          sum(as.genotypes(x)[, nameToKey(x, colnames(m)[j])])
                          )

                    if(has.npb.bootstrap) {
                        df$BOOT.NPB = 
                            c(df$BOOT.NPB,
                              npb.boot.conf[[z]][ rownames(m)[i], colnames(m)[j] ] * 100
                              )
                    } else {
                        df$BOOT.NPB = c(df$BOOT.NPB, NA)
                    }

                    if(has.pb.bootstrap) {
                        df$BOOT.PB = 
                            c(df$BOOT.PB,
                              pb.boot.conf[[z]][ rownames(m)[i], colnames(m)[j] ] * 100
                              )
                    } else {
                        df$BOOT.PB = c(df$BOOT.PB, NA)
                    }

                    if(has.sb.bootstrap) {
                        df$BOOT.SB = 
                            c(df$BOOT.SB,
                              sb.boot.conf[[z]][ rownames(m)[i], colnames(m)[j] ] * 100
                              )
                    } else {
                        df$BOOT.SB = c(df$BOOT.SB, NA)
                    }
                }
            }
        }

        df = cbind(df$SELECTS,
                   df$SELECTED,
                   df$OBS.SELECTS,
                   df$OBS.SELECTED,
                   df$BOOT.NPB,
                   df$BOOT.PB,
                   df$BOOT.SB)

        colnames(df) = 
            c('SELECTS',
              'SELECTED',
              'OBS.SELECTS',
              'OBS.SELECTED',
              'NONPAR.BOOT',
              'PAR.BOOT',
              'STAT.BOOT')

        rownames(df) = paste(1:nrow(df))

        df = data.frame(df, stringsAsFactors = FALSE) 
        df$OBS.SELECTS = as.numeric(df$OBS.SELECTS)
        df$OBS.SELECTED = as.numeric(df$OBS.SELECTED)

        if (!has.npb.bootstrap) {
            df$NONPAR.BOOT = NULL
        } else {
            df$NONPAR.BOOT = as.numeric(df$NONPAR.BOOT)
        }

        if (!has.pb.bootstrap) {
            df$PAR.BOOT = NULL
        } else {
            df$PAR.BOOT = as.numeric(df$PAR.BOOT)
        }

        if (!has.sb.bootstrap) {
            df$STAT.BOOT = NULL
        } else {
            df$STAT.BOOT = as.numeric(df$STAT.BOOT)
        }
        return(df)
    }

    res = lapply(seq_along(matrix), matrix.to.df)
    names(res) = names(matrix)

    return(res)
}


#' Returns a dataframe with all the average/stdev entropy loss score of a 
#' TRONCO model. It is possible to specify models if multiple
#' reconstruction have been performed. 
#'
#' @examples
#' data(test_model)
#' as.kfold.eloss(test_model)
#' as.kfold.eloss(test_model, models='aic')
#'
#' @title as.kfold.eloss
#' @param x A TRONCO model.
#' @param models A subset of reconstructed models, all by default.
#' @param values If you want to see also the values
#' @return All the bootstrap scores in a TRONCO model 
#' @export as.kfold.eloss
#' 
as.kfold.eloss <- function(x,
                           models = names(x$model),
                           values = FALSE) {
    is.compliant(x)
    is.model(x)

    if(is.null(x$kfold)) {
        stop('Cross-validation was never performed on this model!')
    }

    ret = NULL
    ret$C1 = NULL
    ret$C2 = NULL
    ret$C3 = NULL
    ret$C4 = NULL

    for (i in 1:length(models)) {
        if (!is.null(get(models[i], x$kfold)$eloss)) {
            meanll = mean(get(models[i], x$kfold)$eloss)
            ll = get(models[i], x$model)$logLik
            ratio = meanll / abs(ll) * 100

            ret$C1 = c(ret$C1, meanll)
            ret$C2 = c(ret$C2, ratio)
            ret$C3 = c(ret$C3, sd(get(models[i], x$kfold)$eloss))
            ret$C4 = c(ret$C4, paste(round(get(models[i], x$kfold)$eloss, 2), sep = ', ', collapse = ', '))
        }
    }

    ret = data.frame(ret)
    rownames(ret) = models
    colnames(ret) = c('Mean', '%-of-logLik', 'Stdev', 'Values (rounded, 2 digits)')
    
    if(!values) ret = ret[, 1:3]

    return(ret)
}


#' Returns a dataframe with all the prediction error score in a 
#' TRONCO model. It is possible to specify a subset of events
#' or models if multiple reconstruction have been performed. 
#'
#' @examples
#' data(test_model)
#' as.kfold.prederr(test_model)
#' as.kfold.prederr(test_model, models='aic')
#'
#' @title as.kfold.prederr
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param values If you want to see also the values
#' @return All the bootstrap scores in a TRONCO model 
#' @export as.kfold.prederr
#' 
as.kfold.prederr <- function(x,
                             events = as.events(x),
                             models = names(x$model),
                             values = FALSE) {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    if (is.null(x$kfold) ) {
        stop('Crossvalidation was not executed for this object!')
    }

    matrix.to.df <- function(z) {   
        if ( any(is.null(x$kfold[[z]]$prederr)) ) {
            stop('Crossvalidation was not executed for the required model: ', models[z] )
        }

        # Already prepared data - just wrap it in a dataframe
        df = NULL
        df$prederr = t(as.data.frame(x$kfold[[z]]$prederr)) # values

        df$MEAN.PREDERR = apply(df$prederr, 1, mean) # means
        df$SD.PREDERR = apply(df$prederr, 1, sd) # standard deviation
        df$VALUES.PREDERR = apply(df$prederr, 1, function(z){paste(round(z, 3), collapse = ', ')}) # collapse values

        df$SELECTED = gsub("\\.", " ", rownames(df$prederr)) # for later merge we use this

        df = data.frame(df, stringsAsFactors = FALSE) 
        rownames(df) = paste(1:nrow(df))
        df = df[, c('SELECTED', 'MEAN.PREDERR', 'SD.PREDERR', 'VALUES.PREDERR')]
        colnames(df)[4] = 'PREDERR Values (rounded, 3 digits)'

        ## Filter out events if required.

        sel.events = apply(events, 1, function(z){paste(z, collapse = ' ')})
        df = df[which(df$SELECTED %in% sel.events), , drop = F]

        if(!values) df = df[, 1:3]
        
        return(df)
    }

    res = lapply(seq_along(models), matrix.to.df)
    names(res) = models

    return(res)
}


#' Returns a dataframe with all the posterior classification error
#' score in a TRONCO model. It is possible to specify a subset of events
#' or models if multiple reconstruction have been performed. 
#'
#' @examples
#' data(test_model)
#' as.kfold.posterr(test_model)
#' as.kfold.posterr(test_model, events=as.events(test_model)[5:15,])
#'
#' @title as.kfold.posterr
#' @param x A TRONCO model.
#' @param events A subset of events as of \code{as.events(x)}, all by default.
#' @param models A subset of reconstructed models, all by default.
#' @param values If you want to see also the values
#' @return All the posterior classification error scores in a TRONCO model 
#' @export as.kfold.posterr
#' 
as.kfold.posterr <- function(x,
                             events = as.events(x),
                             models = names(x$model),
                             values = FALSE) {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    matrix = 
    as.adj.matrix(x,
        events = events,
        models = models,
        type = 'fit')


    matrix = lapply(matrix, keysToNames, x = x)

    matrix.to.df <- function(z) {  
        if ( is.null(get(z, x$kfold)$posterr) ) {
            return(NULL)
        }

        m = matrix[[z]]
        posterr.matrix = keysToNames(get(z, x$kfold)$posterr, x = x)

        entries = length(which(m == 1))
        df = NULL
        df$SELECTS = NULL
        df$SELECTED = NULL
        df$MEAN.POSTERR = NULL
        df$SD.POSTERR = NULL
        df$POSTERR = NULL

        if (entries == 0) {
            return(NULL)
        }

        for (i in 1:ncol(m)) {
            for (j in 1:nrow(m)) {
                if (m[i, j] == 1) {
                    select = rownames(m)[i]
                    df$SELECTS = c(df$SELECTS, select)
                    selected = colnames(m)[j]
                    val = posterr.matrix[[select,selected]]
                    df$SELECTED = c(df$SELECTED, colnames(m)[j])
                    df$MEAN.POSTERR = c(df$MEAN.POSTERR, mean(val))
                    df$SD.POSTERR = c(df$SD.POSTERR, sd(val))
                    posterr = paste(round(val, 3), collapse=', ')
                    df$POSTERR = c(df$POSTERR, posterr) 
                }
            }
        }

        
        df = cbind(df$SELECTS,
                   df$SELECTED,
                   df$MEAN.POSTERR,
                   df$SD.POSTERR,
                   df$POSTERR)

        rownames(df) = c(1:nrow(df))
        df = data.frame(df, stringsAsFactors = FALSE) 

        colnames(df) = 
            c('SELECTS',
              'SELECTED',
              'MEAN.POSTERR',
              'SD.POSTERR',
              'POSTERR Values (rounded, 3 digits)')

        if(!values) df = df[, 1:4]
        

        return(df)
    }

    res = lapply(models, matrix.to.df)
    names(res) = names(matrix)

    return(res)
}


# Get parents for each node
#
# @title as.parents.pos
# @param x A TRONCO model.
# @param events A subset of events as of \code{as.events(x)}, all by default.
# @param models A subset of reconstructed models, all by default.
# @return A list of parents for each node
#
as.parents.pos <- function(x,
                           events = as.events(x),
                           models = names(x$model)) {
    is.compliant(x)
    is.model(x)
    is.events.list(x, events)

    m = as.models(x, models = models)

    ret = list()
    for (i in models) {
        mat = m[[i]]$parents.pos
        mat = mat[rownames(events), , drop = FALSE]
        ret = append(ret, list(mat)) 
    }

    names(ret) = models
    return(ret) 
}


#' Get parameters of a model
#'
#' @title as.parameters
#' @param x A TRONCO model.
#' @return A list of parameters
#' @export as.parameters
#
as.parameters <- function(x) {
    is.compliant(x)
    is.model(x)

    return(x$parameters) 
}


#' Remove models, hypothesis and patterns from a TRONCO model.
#'
#' @title as.dataset
#' @param x A TRONCO model.
#' @return A TRONCO model without models, hypothesis and patterns.
#' @export as.dataset
#'
as.dataset <- function(x) {
    is.compliant(x)
    if (has.model(x)) {
        ret = delete.model(x)
    }
    for (pattern in as.patterns(ret)) {
        ret = delete.pattern(ret, pattern)
    }

    return(ret) 
}


#' Return true if the TRONCO dataset 'x', which should be a TRONCO compliant dataset 
#' - see \code{is.compliant} - has stage annotations for samples. Some sample stages 
#' might be annotated as NA, but not all.
#'
#' @examples
#' data(test_dataset)
#' has.stages(test_dataset)
#' data(stage)
#' test_dataset = annotate.stages(test_dataset, stage)
#' has.stages(test_dataset)
#'
#' @title has stages
#' @param x A TRONCO compliant dataset.
#' @return TRUE if the TRONCO dataset has stage annotations for samples.
#' @export has.stages
#' 
has.stages <- function(x) {
    return(!(all(is.null(x$stages)) || all(is.na(x$stages))))
}


#' Return true if there are duplicated events in the TRONCO dataset 'x', which should be
#' a TRONCO compliant dataset - see \code{is.compliant}. Events are identified by a gene 
#' name, e.g., a HuGO_Symbol, and a type label, e.g., c('SNP', 'KRAS')
#'
#' @examples
#' data(test_dataset)
#' has.duplicates(test_dataset)
#' 
#' @title has.duplicates
#' @param x A TRONCO compliant dataset.
#' @return TRUE if there are duplicated events in \code{x}.
#' @export has.duplicates
#' 
has.duplicates <- function(x) {  
    ## Find duplicate over the dataset.
    dup = duplicated(as.events(x))

    ## Return true if at least one duplicate is found.
    return(any(dup))
}


#' Return true if there is a reconstructed model in the TRONCO dataset 'x', which should be
#' a TRONCO compliant dataset - see \code{is.compliant}.
#'
#' @examples
#' data(test_dataset)
#' has.model(test_dataset)
#' 
#' @title has.model
#' @param x A TRONCO compliant dataset.
#' @return TRUE if there is a reconstructed model in \code{x}.
#' @export has.model
#' 
has.model <- function(x) {
    is.compliant(x)
    if (length(x$model) > 0)            # Isn't 'return(length(x$mode) > 0)' enough?
        return(TRUE)
    return(FALSE)
}


#' Return the events duplicated in \code{x}, if any. Input 'x' should be
#' a TRONCO compliant dataset - see \code{is.compliant}. 
#'
#' @examples
#' data(test_dataset)
#' duplicates(test_dataset)
#'
#' @title duplicates
#' @param x A TRONCO compliant dataset.
#' @return A subset of \code{as.events(x)} with duplicated events.
#' @export duplicates
#' 
duplicates <- function(x) {
    is.compliant(x)
    return(as.events(x)[duplicated(as.events(x)),])
}


#' Print to console a short report of a dataset 'x', which should be
#' a TRONCO compliant dataset - see \code{is.compliant}. 
#'
#' @examples
#' data(test_dataset)
#' view(test_dataset)
#' 
#' @title view
#' @param x A TRONCO compliant dataset.
#' @param view The firse \code{view} events are shown via \code{head}.
#' @export view
#' 
view <- function(x, view = 5) {
    is.compliant(x)
    x = enforce.numeric(x)
    view = min(view, nevents(x))

    if (as.description(x) != "")
        cat(paste('Description: ', as.description(x), '.\n', sep = ''))

    cat(paste('-- TRONCO Dataset: n=',
              nsamples(x),
              ', m=',
              nevents(x),
              ', |G|=',
              ngenes(x),
              ', patterns=',
              npatterns(x),
              '.\n',
              sep = '')) 
    cat(paste('Events (types): ',
              paste(as.types(x), collapse=', '),
              '.\n',
              sep = ''))
    cat(paste('Colors (plot): ',
              paste(as.colors(x), collapse=', '),
              '.\n',
              sep = ''))

    if (has.stages(x)) {
        cat(paste('Stages: '))

        s = unlist(sort(unique(as.stages(x)[, 1])))
        cat((paste(paste(s, collapse=', ', sep = ''), '.\n', sep = '')))
    }

    cat(paste('Events (', view, ' shown):\n', sep = ''))
    to.show =
        paste( '\t',
              rownames(as.events(x)[1: view,]), ':',
              as.events(x)[1: view, 1],
              as.events(x)[1: view, 2],
              sep = ' ')

    cat(paste(to.show, collapse = '\n'))

    cat(paste('\nGenotypes (', view, ' shown):\n', sep = ''))
    print(head(x$genotypes[,1:view, drop = FALSE]))

    if (has.model(x)) {
        cat('\n-- TRONCO Model(s): ', x$parameters$algorithm, '\n')

        if (x$parameters$algorithm == 'CAPRI') {
            cat('Score optimization via ')
            if (x$parameters$command == 'hc') {
                cat('Hill-Climbing with ')
            }

            if (x$parameters$command == 'tabu') {
                cat('Tabu Search with ')
            }
            
            cat(paste(toupper(x$parameters$regularization), sep = ', ', collapse = ', '),' regularizers.\n', sep = '')
        
            if ('bic' %in% x$parameters$regularization) {
                cat('BIC: ')
                cat('score', x$model$bic$score ,'|',
                    'logLik', x$model$bic$logLik ,'|', 
                    nrow(as.selective.advantage.relations(x)$bic), 'selective advantage relations.\n')
            }
          
            if ('aic' %in% x$parameters$regularization) {
                cat('AIC: ')
                cat('score', x$model$aic$score ,'|',
                    'logLik', x$model$aic$logLik ,'|', 
                    nrow(as.selective.advantage.relations(x)$aic), 'selective advantage relations.\n')
            }
        }
      
        cat('Available confidence measures:\n')
        cat('\tp-values: Temporal priority | Probability raising | Hypergeometric\n')

        if (!(is.null(x$confidence) && is.na(x$confidence))) {
            models = names(x$model)
        }
      
        has.npb.bootstrap = is.null(x$bootstrap[[models[1]]]$npb)
        has.pb.bootstrap = is.null(x$bootstrap[[models[1]]]$pb)
        has.sb.bootstrap = is.null(x$bootstrap[[models[1]]]$sb)
      
        if(!is.null(x$bootstrap)) {
            cat('\tBootstrap scores: ')
            if (!is.null(x$bootstrap[[models[1]]]$npb)) {
                cat('Non-parametric | ')
            }
            if (!is.null(x$bootstrap[[models[1]]]$pb)) {
                cat('Parametric | ')
            }
            if (!is.null(x$bootstrap[[models[1]]]$sb)) {
                cat('Statistical')
            }
        }
        cat('\n')
    }
} 



#' Return the number of types in the dataset.
#'
#' @examples
#' data(test_dataset)
#' ntypes(test_dataset)
#'
#' @title ntypes
#' @param x A TRONCO compliant dataset.
#' @return The number of types in the dataset.
#' @export ntypes
#' 
ntypes <- function(x) {
    return(length(as.types(x)))
}


#' Return the number of samples in the dataset.
#'
#' @examples
#' data(test_dataset)
#' nsamples(test_dataset)
#'
#' @title nsamples
#' @param x A TRONCO compliant dataset.
#' @return The number of samples in the dataset.
#' @export nsamples
#' 
nsamples <- function(x) {
    return(nrow(x$genotypes))
}


#' Return the number of events in the dataset involving a certain gene or type of event.
#'
#' @examples
#' data(test_dataset)
#' nevents(test_dataset)
#'
#' @title nevents
#' @param x A TRONCO compliant dataset.
#' @param genes The genes to consider, if NA all available genes are used.
#' @param types The types of events to consider, if NA all available types are used.
#' @return The number of events in the dataset involving a certain gene or type of event.
#' @export nevents
#' 
nevents <- function(x, genes = NA, types = NA) {
    return(nrow(as.events(x, genes, types)))
}


#' Return the number of genes in the dataset involving a certain type of event.
#'
#' @examples
#' data(test_dataset)
#' ngenes(test_dataset)
#'
#' @title ngenes
#' @param x A TRONCO compliant dataset.
#' @param types The types of events to consider, if NA all available types are used.
#' @return The number of genes in the dataset involving a certain type of event.
#' @export ngenes
#' 
ngenes <- function(x, types = NA) {
    return(length(as.genes(x, types = types)))
}


#' Return the number of patterns in the dataset
#'
#' @examples
#' data(test_dataset)
#' npatterns(test_dataset)
#'
#'
#' @param x the dataset.
#' @export npatterns
#' 
npatterns <- function(x) {
    if (any(is.null(x$hypotheses))) {
        return(0)
    }

    if ('hstructure' %in% names(x$hypotheses)) {
        return(length(ls(x$hypotheses$hstructure)))
    }
    return(0)
}


#' Return the number of hypotheses in the dataset
#'
#' @examples
#' data(test_dataset)
#' nhypotheses(test_dataset)
#'
#' @param x the dataset.
#' @export nhypotheses
#' 
nhypotheses <- function(x) {
    if (npatterns(x) < 1) {
        return(0)
    }

    if ('hlist' %in% names(x$hypotheses)) {
        return(length(x$hypotheses$hlist) / 2)
    }
    return(0)
}


#' Convert the internal reprensentation of genotypes to numeric, if not. 
#'
#' @examples
#' data(test_dataset)
#' test_dataset = enforce.numeric(test_dataset)
#'
#' @title enforce.numeric
#' @param x A TRONCO compliant dataset.
#' @return Convert the internal reprensentation of genotypes to numeric, if not.
#' @export enforce.numeric
#' 
enforce.numeric <- function(x) {
    is.compliant(x)
    if (!all(is.numeric(x$genotypes[1, ]))) {
        rn = as.samples(x)
        x$genotypes = apply(x$genotypes, 2, as.numeric)
        rownames(x$genotypes) = rn
    }

    return(x)
}


#' Convert the internal representation of genotypes to character, if not. 
#'
#' @examples
#' data(test_dataset)
#' test_dataset = enforce.string(test_dataset)
#'
#' @title enforce.string
#' @param x A TRONCO compliant dataset.
#' @return Convert the internal reprensentation of genotypes to character, if not.
#' @export enforce.string
#' 
enforce.string <- function(x) {
    is.compliant(x)
    if (!all(is.character(x$genotypes[1,]))) {
        rn = as.samples(x)
        x$genotypes = apply(x$genotypes, 2, as.character)
        rownames(x$genotypes) = rn
    }

    return(x)
}


#' Sort the internal genotypes according to event frequency.
#'
#' @examples
#' data(test_dataset)
#' sort.by.frequency(test_dataset)
#'
#' @title sort.by.frequency
#' @param x A TRONCO compliant dataset.
#' @param decreasing Inverse order. Default TRUE
#' @param ... just for compatibility
#' @return A TRONCO compliant dataset with the internal genotypes sorted according to event frequency.
#' @export sort.by.frequency
#' 
sort.by.frequency <- function(x, decreasing = TRUE, ...) {
    other.argument = list(...)

    is.compliant(x)

    x = enforce.numeric(x)
    sums = colSums(x$genotypes)

    x$genotypes = x$genotypes[, order(sums, decreasing = decreasing), drop = F]
    x$annotations = x$annotations[colnames(x$genotypes), , drop = F]

    return(x)  
}


#' Convert colnames/rownames of a matrix into intelligible event names, e.g., change a key G23 in 'Mutation KRAS'.
#' If a name is not found, the original name is left unchanged.
#'
#' @examples
#' data(test_model)
#' adj_matrix = as.adj.matrix(test_model, events=as.events(test_model)[5:15,])$bic
#' keysToNames(test_model, adj_matrix)
#'
#'
#' @title keysToNames
#' @param x A TRONCO compliant dataset.
#' @param matrix A matrix with colnames/rownames which represent genotypes keys. 
#' @return The matrix with intelligible colnames/rownames. 
#' @export keysToNames
#' 
keysToNames <- function(x, matrix) {
    is.compliant(x)
    if (!is.matrix(matrix)
        || any(is.null(colnames(matrix)))
        || any(is.null(rownames(matrix))))
        stop('"matrix" should be a matrix with rownames/colnames.')

    events = as.events(x)
    resolve <- function(y) { 
        if (y %in% rownames(events)) paste(events[y,], collapse=' ')
        else y
    }

    colnames(matrix) = sapply(colnames(matrix), resolve)
    rownames(matrix) = sapply(rownames(matrix), resolve)
    return(matrix)
}


nameToKey <- function(x, name) {
    is.compliant(x)

    types = as.types(x)
    for (i in types) {
        if (nchar(name) > nchar(i) && substr(name, 1, nchar(i)) == i)
            return(rownames(as.events(x, 
                                      genes = substr(name, nchar(i) + 2, nchar(name)),
                                      types = substr(name, 1, nchar(i)))))
    }
    
    stop('"name" is not a key!')
}


# Check if logic node down
#
# @title is logical node down
# @param node A node identifier
# @return boolean
#
is.logic.node.down <- function(node) {
    if (substr(node, start = 1, stop = 3) == 'OR_')
        return(TRUE)
    if (substr(node, start = 1, stop = 4) == 'XOR_')
        return(TRUE)
    if (substr(node, start = 1, stop = 4) == 'AND_')
        return(TRUE)
    if (substr(node, start = 1, stop = 4) == 'NOT_')
        return(TRUE)
    return(FALSE)
}


# Check if logic node up
#
# @title is logical node up
# @param node A node identifier
# @return boolean
is.logic.node.up <- function(node) {
    if (substr(node, start = 1, stop = 2) == 'UP')
        return(TRUE)
    return(FALSE)
}


# Check if logic node down or up
#
# @title is logical node
# @param node A node identifier
# @return boolean

is.logic.node <- function(node) {
    return(is.logic.node.up(node) || is.logic.node.down(node))
}


#### end of file -- as.functions.R
