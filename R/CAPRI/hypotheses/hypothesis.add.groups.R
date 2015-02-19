# commenti

hypothesis.add.group = function(x, FUN, group, ...) {
  op = deparse(substitute(FUN))
  
  effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
  effect = paste(effect, collapse = ', ')
  
  ngroup = length(group)
  if(ngroup < 2) 
    return()
  
  hom.group = lapply(group, function(g, x){if(nevents(x, genes = g) > 1) T else F }, x)
  hom.group = group[unlist(hom.group)]
  
  gene.hom = function(g, h) {
    if (g %in% h)
      return(paste0('OR(\'', g, '\')'))
    return(paste0('\'', g, '\''))
  }
  
  # print(hom.group)
  
  tot = 0
  for(k in 2:ngroup) { tot = tot + (factorial(ngroup) / factorial(k))}
  
  cat('*** Number of hypothses to be generated: ', tot, '\n')
  if (length(hom.group) > 0)
    cat('*** Genes with functional homologous found: ', unlist(hom.group), '\n')
  
  # create a progress bar
  pb <- txtProgressBar(1, tot, style = 3)
  pbPos = 1
  
  e = list()
    
  for(i in 2:ngroup) {
    gr = combn(unlist(group), i)
    
    for(j in 1:ncol(gr))  {	
      genes = as.list(gr[,j])

      ### serve per i singoli con + mutazioni
      #if (length(genes) == 1 && nevents(x, genes = genes[1]) < 2) {
      #  next
      #}
      
      #start the progress bar
      setTxtProgressBar(pb, pbPos)
      pbPos = pbPos + 1
        
      hypo.name = paste(unlist(genes), sep='_', collapse='_')
      hypo.genes = paste(lapply(genes, function(g, hom.group){gene.hom(g, hom.group)}, hom.group), collapse=', ')  
      
      hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                        op, '_', hypo.name, 
                        '\', lifted.formula = ',
                        op,
                        '(',
                        hypo.genes,
                        '), ',
                        effect,
                        ')')
      
      # cat('*** Evaluating ', hypo.add, '\n')
      
      tryCatch({
        x = eval(parse(text=hypo.add))
      },
      error=function(cond){
        m = paste('Error on', hypo.add, '.\n', cond)
        e = append(e, m)
        message(e)
      },
      warning=function(cond){
        m = paste('Warning on', hypo.add, '.\n', cond)
        e = append(e, m)
        message(e)
      })
      
    }
    
  }
  
  # close progress bar
  close(pb)
  
  return(x)
}


hypothesis.add.hom = function(x, ..., genes = as.genes(x)){
  # in questa funzione, per ogni gene che ha più di un tipo di alterazione
  # aggiungo l'OR
  
  effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
  effect = paste(effect, collapse = ', ')
  
  hom.group = lapply(genes, function(g, x){if(nevents(x, genes = g) > 1) T else F }, x)
  hom.group = genes[unlist(hom.group)]
  
  cat('*** Genes with functional homologous found: ', hom.group, '\n')
  
  # create a progress bar
  pb <- txtProgressBar(1, length(hom.group), style = 3);
  
  for(i in 1:length(hom.group)) {
    
    #start the progress bar
    setTxtProgressBar(pb, i)
    
    hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                      'OR_', hom.group[[i]], 
                      '\', lifted.formula = OR(\'',
                      hom.group[[i]],
                      '\'), ',
                      effect,
                      ')')
          
    tryCatch({
      x = eval(parse(text=hypo.add))
    },
    error=function(cond){
      message(paste('Error on', hypo.add, '.'))
      message(cond)
      message('\n')
    },
    warning=function(cond){
      message(paste('Warning on', hypo.add, '.'))
      message(cond)
    })
  }
  
  # close progress bar
  close(pb)
  
  return(x)
  
}

