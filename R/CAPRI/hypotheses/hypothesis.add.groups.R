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
      print(genes)
      formula = quote(FUN(unlist(genes)))
      print(formula)
      
      x = hypothesis.add(x, label.formula = 'asd', lifted.formula = formula, '*')
      
    }
    cat('\n\n')
    
  }
  
  return()
}