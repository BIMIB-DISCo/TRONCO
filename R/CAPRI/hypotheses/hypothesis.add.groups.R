# commenti

hypothesis.add.groups = function(x, FUN, group, ...) {
  op = deparse(substitute(FUN))
  
  effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
  effect = paste(effect, collapse = ', ')
  
  ngroup = length(group)
  if(ngroup < 2) 
    return()
  
  for(i in 2:ngroup) {
    gr = combn(unlist(group), i)
    
    for(j in 1:ncol(gr))  {	
      genes = as.list(gr[,j])

      ### serve per i singoli con + mutazioni
      #if (length(genes) == 1 && nevents(x, genes = genes[1]) < 2) {
      #  next
      #}
        
      hypo.name = paste(unlist(genes), sep='_', collapse='_')
      hypo.genes = paste(unlist(genes), collapse='\', \'')  
      
      hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                        op, '_', hypo.name, 
                        '\', lifted.formula = ',
                        op,
                        '(\'',
                        hypo.genes,
                        '\'), ',
                        effect,
                        ')')
      
      cat('*** Evaluating ', hypo.add, '\n')
      
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
    
  }
  
  return(x)
}


hypothesis.add.hom = function(x, ..., genes = as.genes(x)){
  # in questa funzione, per ogni gene che ha più di un tipo di alterazione
  # aggiungo l'OR
  
  effect = sapply(as.list(substitute(list(...)))[-1L], deparse)
  effect = paste(effect, collapse = ', ')
  
  # create a progress bar
  pb <- txtProgressBar(1, length(genes), style = 3);
  
  for(i in 1:length(genes)) {
    
    #start the progress bar
    setTxtProgressBar(pb, i)
    
    if(nevents(x, genes = genes[[i]]) > 1) {
      hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                        'OR_', genes[[i]], 
                        '\', lifted.formula = OR(\'',
                        genes[[i]],
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
  }
  
  # close progress bar
  close(pb);
  
  return(x)
  
}

