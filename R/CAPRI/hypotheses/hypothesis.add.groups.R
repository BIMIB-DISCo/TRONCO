# commenti

hypothesis.add.groups = function(x, FUN, group, ...) {
  op = deparse(substitute(FUN))
  
  effect = list(...)

  ngroup = length(group)
  if(ngroup < 2) 
    return()
  
  for(i in 2:ngroup) {
    cat('i', i)
    gr = combn(unlist(group), i)
    
    for(j in 1:ncol(gr))
    {	

      cat('\n')
      genes = as.list(gr[,j])
      
      #formula = quote(FUN(unlist(genes)))
      #print(formula)
      #formula = substitute(f(arg), list(f = FUN, arg = unlist(genes)))
      #print(formula)
      
      hypo.name = paste(unlist(genes), sep='_', collapse='_')
      hypo.genes = paste(unlist(genes), collapse='\', \'')
      hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                        op, '_', hypo.name, 
                        '\', lifted.formula = ',
                        op,
                        '(\'',
                        hypo.genes,
                        '\'), \'*\' )')
      
      print('hadd')
      print(hypo.add)
      
      tryCatch({
        x = eval(parse(text=hypo.add))
      },
      error=function(cond){
        message(paste('Error on', hypo.add, '.\n'))
        message(cond)
        message('\n')
      },
      warning=function(cond){
        essage(paste('Warning on', hypo.add, '.\n'))
        message(cond)
        essage('\n')
      })
      
      
      #eval('hypothesis.add(s, label.formula = \'asd\', lifted.formula = XOR(\'APC\', \'KRAS\'), c(\'CUBN\', \'SNV\'))')
      
      #x = hypothesis.add(x, label.formula = 'asd', lifted.formula = formula, c('CUBN', 'SNV'))
      
    }
    cat('\n\n')
    
  }
  
  return(x)
}