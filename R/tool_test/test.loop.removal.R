#the local path of my git directory
my.GIT = '/Users/daniele/Documents/BIMIB/'

#set the working directory
work.dir = paste0(my.GIT,'TRONCO/R/tool_test');
setwd(work.dir);

#set git directory
invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/CAPRI'),full.names=TRUE,recursive=TRUE),source));
invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/CAPRESE'),full.names=TRUE,recursive=TRUE),source));
invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/TRONCO'),full.names=TRUE,recursive=TRUE),source));
invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/utilities'),full.names=TRUE,recursive=TRUE),source));

#file where the dataset is saved
file.dataset = paste0(work.dir,'/data/test.loop.removal.Rdata');

#load and save the dataset
load(file.dataset);
data = delete.gene(lift, "ACVR1B");
data$hypotheses$patterns = data$hypotheses$patterns[-which(names(data$hypotheses$patterns) == "XOR_ACVR1B")];

#run CAPRI with bic
capri.with.bootstrap.bic = tronco.capri(data);

#plot the results
tronco.plot(capri.with.bootstrap.bic);

#run CAPRI with aic
capri.with.bootstrap.aic = tronco.capri(data,REGULARIZATION="aic");

#plot the results
tronco.plot(capri.with.bootstrap.aic);
