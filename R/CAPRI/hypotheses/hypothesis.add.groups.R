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
    
    for(j in 1:ncol(gr))
    {	
      genes = as.list(gr[,j])
        
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