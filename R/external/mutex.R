######## Export to mutex

## TODO


######## Import mutex groups -- current Mutex version is XXX
# Create a list of unique Mutex groups for a fdr cutoff
# (ref: http://.....)
#' @param file Mutex results ("ranked-groups.txt" file)
#' @param fdr cutoff for fdr
#' @param display print summary table of extracted groups
import.mutex.groups = function(file, fdr=.2, display = TRUE)
{
  # Found somewhere on the web - makes sense
  read.irregular <- function(filenm) 
  {
    fileID <- file(filenm,open="rt")
    nFields <- count.fields(fileID)
    mat <- matrix(nrow=length(nFields),ncol=max(nFields))
    invisible(seek(fileID,where=0,origin="start",rw="read"))
    for(i in 1:nrow(mat) ) {
      mat[i,1:nFields[i]] <-scan(fileID,what="",nlines=1,quiet=TRUE)
    }
    close(fileID)
    df <- data.frame(mat, stringsAsFactors=FALSE)
    return(df)
  }
  
  x = read.irregular(file)
  
  # Check header
  if(any(x[1,1:3] != c('Score', 'q-val', 'Members')))
    warning('File header does not seem to contain \'Score\', \'q-val\' and \'Members field\' - are you
            sure this is a Mutex result file?' )

  # Remove header
  cat(paste('*** Groups extracted - ', (nrow(x) -1), ' total groups.\n', sep=''))
  x = x[-1, , drop = F] # this is c('Score', 'q-val', 'Members')
  x[, 1] = as.numeric(x[,1]) # fdr
  x[, 2] = as.numeric(x[,2]) # q-value
  
  # remove groups  with low fdr
  res = x[which(x[,1] < fdr), , drop = F] 
  
  # remove duplicated groups (permutations)
  res.g = res[, 3:ncol(res)]
  
  for(i in 1:nrow(res.g)) res[i,3:ncol(res)] = sort(res.g[i,], na.last = T)   
    res = res[!duplicated((res[,3:ncol(res), drop=FALSE])), ] 
  
  cat(paste('Selected ', nrow(res), ' unique groups with fdr < ', fdr, '\n', sep=''))
  
  # Create groups
  groups = function(g) {
  	# print(g)
  	
  	# print(g[3:length(g)])
  	
    g = g[3:length(g)]
    g = g[!is.na(g)]
    names(g) = NULL
    
    # print(g)
    return(sort(g))
  }
  
  # print(res)
  # print(apply(res, 1, groups))
  
  # G = apply(res, 1, groups)
  # print(G)
  
  # print(res)
  G = list()
  for(i in 1:nrow(res))
  {
  	gr = list(groups(res[i, ]))
  	names(gr) = paste('MUTEX_GROUP', i, sep='')
  	G = append(G,gr)	
  }
  # print(G)
  
  # if(!is.list(G)) 
  # {
  	# G = list(as.vector(G))	
    # names(G) = paste('MUTEX_GROUP', 1, sep='')
  # }

  # names(G) = paste('MUTEX_GROUP', 1:length(names(G)), sep='')
  rownames(res) = names(G)
  colnames(res)[1:2] = c('fdr', 'score')
  
  # Summary report
  if(display) 
  {	
  	print(res)
 	# print(G)
 	}
  return(G)
}
