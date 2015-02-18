# commenti

hypothesis.add.groups = function(x, FUN, ...) {
  op = deparse(substitute(FUN))
  
  group = list(...)

  ngroup = length(group)

  if(ngroup < 2) 
    return()
  
  for(i in 2:ngroup) {
    cat('i', i)
    gr = combn(unlist(group), i)
    
    for(j in 1:ncol(gr))
    {	
      cat('\n')
      
      genes = list()
      
      for(k in 1:nrow(gr)) {
        cat(gr[k,j], '')
        genes = append(genes, gr[k,j])
      }
      cat('\n')
      #print(genes)
      formula = quote(FUN(unlist(genes)))
      #print(formula)
      #formula = substitute(f(arg), list(f = FUN, arg = unlist(genes)))
      print(formula)
      
      hypo.name = paste(op, unlist(genes), sep='_')
      hypo.genes = paste(unlist(genes), sep='\', \'')
      hypo.add = paste0('hypothesis.add(x, label.formula = \'', 
                        hypo.name, 
                        '\', lifted.formula = ',
                        op,
                        '(\'',
                        hypo.genes,
                        '\'), \'*\' )')
      
      print('hadd')
      print(hypo.add)
      eval(parse(text=hypo.add))
      
      #eval('hypothesis.add(s, label.formula = \'asd\', lifted.formula = XOR(\'APC\', \'KRAS\'), c(\'CUBN\', \'SNV\'))')
      
      #x = hypothesis.add(x, label.formula = 'asd', lifted.formula = formula, c('CUBN', 'SNV'))
      
    }
    cat('\n\n')
    
  }
  
  return()
}