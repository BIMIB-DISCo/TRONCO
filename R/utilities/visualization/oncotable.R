# PDF and laex tables
genes.table.report = function(x, name, dir=getwd(), maxrow=33, 
                        font=10, height=11, width=8.5, fill="lightblue") 
{
  # Print table with gridExtra and xtables
  print.table = function(table, name, dir=getwd(), maxrow, font, height, 
                         width, fill)
  {
    cat('Printing PDF and Latex table to files: \n')
    cat(paste('PDF \t\t', dir, '/', name, '.genes-table.pdf\n', sep=''))
    cat(paste('Latex\t\t', dir, '/', name, '.genes-table.tex\n', sep=''))
    
    # output pdf
    require(gridExtra)  
    require(xtable)  
    
    cur.dev = dev.cur()
    
    pdf(file=paste(dir, '/', name, '.genes-table.pdf', sep=''), height=height, width=width)
    
    # max rows per page	
    npages = ceiling(nrow(table)/maxrow); 
    
    pb = txtProgressBar(1, npages, style = 3);      
    for (i in 1:npages) 
    {
      setTxtProgressBar(pb, i)  
      idx = seq(1+((i-1)*maxrow), i*maxrow); 
      
      if(max(idx) > nrow(table)) idx = idx[idx < nrow(table)]		
      
      grid.newpage(); 
      grid.table(table[idx, ],
                 gpar.coretext = gpar(fontsize = font),
                 gpar.corefill = gpar(fill = fill, alpha=0.5, col = NA),
                 h.even.alpha = 0.5)
    } 
    close(pb)
    
    # output latex    
    print(xtable(table, digits=0), file=paste(dir, '/', name, '.genes-table.tex', sep=''), type='latex')
  
    dev.off()
    dev.set(which=cur.dev)
    }
  
  cat(paste('Preparing output table with ', ngenes(x),' genes ...\n'))
  genes = as.genes(x)
  types = as.types(x)
  
  data = matrix(0, nrow=ngenes(x), ncol=ntypes(x))

  genes.table = data.frame(data, row.names=genes, stringsAsFactors=FALSE)
	colnames(genes.table) = types

  x = enforce.numeric(x)
  
  pb = txtProgressBar(1, ngenes(x), style = 3);
	for(i in 1:ngenes(x))
	{
	  setTxtProgressBar(pb, i)  
		g = as.gene(x, genes=genes[i])
        
    if(ncol(g) > 0)
		{
			gg = colSums(apply(g, 2, as.numeric))
	
			genes.table[rownames(genes.table)[i], colnames(g)] = gg
			genes.table[rownames(genes.table)[i], 'Alterations'] = paste( round(sum(gg) / nsamples(x) * 100), '%', sep='')	
			genes.table[rownames(genes.table)[i], 'Frequency'] = sum(gg) / nsamples(x) 	
		}

	}
  # close progress bar
  close(pb)

  
	genes.table = genes.table[order(genes.table$Frequency, decreasing = TRUE), ]
	genes.table$Frequency = NULL

	print.table(table=genes.table, name=name, dir=getwd(), maxrow=maxrow, font=font, height=height, 
  width=width, fill=fill)
  
	return(genes.table)
}



# stacked histogram
genes.table.plot = function(x, name, dir=getwd()) 
{  
  
  require(reshape2)
  require(ggplot2)
  
  cat('Preparing output table: creating alterations profiles and selecting events with minimum frequency.\n') 
  alterations = sort.by.frequency(as.alterations(x))
  
  cat('Stacked histogram with genes in the following order (head): ')
  cat(head(as.genes(alterations)))
  cat('\nRetrieving all events for the above genes.\n')

  y = events.selection(x, filter.in.names=as.genes(alterations))
  y = enforce.numeric(y)
#  show(y)
  
 # print(ntypes(y))
#  print(ngenes(y))

  table = matrix(rep(0, ntypes(y) * ngenes(y)), ncol=ntypes(y), nrow=ngenes(y))
  
  rownames(table) = 1:ngenes(alterations)
  colnames(table) = as.types(y)

  cat('Populating histogram table.\n')
  pb = txtProgressBar(1, ngenes(alterations), style = 3);      
  for(i in 1:ngenes(alterations))
  {
    setTxtProgressBar(pb, i)  
    
    g = colSums(as.gene(y, gene=as.genes(alterations)[i]))
    table[i, ] = g
  }
  close(pb)

  table = cbind(Rank=as.genes(alterations), table)
  table.melt = melt(as.data.frame(table), id.var='Rank')  

  #print(table)
  #print(table.melt)
  

# Problem, does not work well - can't assign colors in as.colors(y)
  p = ggplot(table.melt, aes(x = Rank, y = value, fill = variable)) 
  p + geom_bar(stat = "identity")
  print(p)

  return(table.melt)
}
# # cat('Latex tables (genes)')
# genes.table.input = genes.report(input)
# genes.table.sub1 = genes.report(sub1)
# genes.table.sub2 = genes.report(sub2)
# genes.table.sub3 = genes.report(sub3)
# genes.table.sub4 = genes.report(sub4)
# genes.table.nh = genes.report(non.hyper)
# 
# idx = rownames(genes.table.input)
# full.table = cbind(TOT=genes.table.input[idx, 'Freq'], genes.table.sub1[idx, ])
# full.table = cbind(
# full.table, genes.table.sub2[idx, c('Amps', 'Dels', 'SNVs', 'Freq')], 
# genes.table.sub3[idx, c('Amps', 'Dels', 'SNVs', 'Freq')], 
# genes.table.sub4[idx, c('Amps', 'Dels', 'SNVs','Freq')])
# 
# print.table(full.table, 'Full-table', font = 10, width = 20)
# print.table(pathway.genes.df, 'genes-list')


