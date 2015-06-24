#the local path of my git directory
my.GIT = '/Users/daniele/Documents/BIMIB/'

#set the working directory
work.dir = paste0(my.GIT,'TRONCO/R/tool_test');
setwd(work.dir);

#file where the dataset is saved
file.dataset.ovarian = paste0(work.dir,'/data/ovarian-data.txt');

#set git directory
# invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/CAPRI'),full.names=TRUE,recursive=TRUE),source));
# invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/CAPRESE'),full.names=TRUE,recursive=TRUE),source));
# invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/TRONCO'),full.names=TRUE,recursive=TRUE),source));
# invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R/utilities'),full.names=TRUE,recursive=TRUE),source));
invisible(sapply(list.files(pattern="[.]R$",path=paste0(my.GIT,'TRONCO/R'),full.names=TRUE,recursive=FALSE),source));

#load the dataset and set all the values for colnames and rownames
genotypes = read.table(file.dataset.ovarian);
colnames(genotypes) = gsub("V", "gene ", colnames(genotypes));
rownames(genotypes) = paste("patient ",rownames(genotypes),sep="");
annotations = array("",c(ncol(genotypes),2));
colnames(annotations) = c("type","event");
rownames(annotations) = colnames(genotypes);
annotations[c(1,2,6),1] = "Gain";
annotations[c(3,4,5,7),1] = "Loss";
annotations[1,2] = "8q+";
annotations[2,2] = "3q+";
annotations[3,2] = "5q-";
annotations[4,2] = "4q-";
annotations[5,2] = "8p-";
annotations[6,2] = "1q+";
annotations[7,2] = "Xp-";
types = array("",c(length(unique(annotations[,"type"])),1));
colnames(types) = "color";
rownames(types) = unique(annotations[,"type"]);
types[1,1] = "Red";
types[2,1] = "Blue";

#create the input variable
data = list(genotypes = genotypes, annotations = annotations, types = types);

#perform the reconstruction with CAPRESE without estimations
caprese.no.estimations = tronco.caprese(data);

#perform the estimations given the error rates
caprese.no.estimations.with.error.rates = tronco.estimation(caprese.no.estimations,list(error.fp=0,error.fn=0));

#perform the estimations without giving the error rates
caprese.no.estimations.without.error.rates = tronco.estimation(caprese.no.estimations);

#perform the reconstruction with CAPRESE with estimations
caprese = tronco.caprese(data,do.estimation=TRUE);

#perform the estimation by non-parametric and parametric bootstraps with CAPRESE
set.seed("12345");
caprese.non.parametric = tronco.bootstrap(caprese);
set.seed("12345");
caprese.parametric = tronco.bootstrap(caprese,type="parametric");

#perform the reconstruction with CAPRI
data = hypothesis.add(data,"H1",OR(XOR(c("8q+","Gain"),c("4q-","Loss")),AND(c("3q+","Gain"),c("5q-","Loss")),c("8p-","Loss"),c("1q+","Gain")),hypothesis.lifted.effects(c("Xp-","Loss")));
data = hypothesis.add(data,"H2",AND(XOR(c("8q+","Gain"),c("4q-","Loss")),OR(c("3q+","Gain"),c("5q-","Loss")),c("8p-","Loss"),c("1q+","Gain")),hypothesis.lifted.effects(c("Xp-","Loss")));
data = hypothesis.add(data,"H3",OR(XOR(c("8q+","Gain"),c("4q-","Loss")),c("8p-","Loss")),hypothesis.lifted.effects("*"));
data = hypothesis.add(data,"H4",OR(AND(c("8q+","Gain"),c("4q-","Loss")),c("8p-","Loss")),hypothesis.lifted.effects(c("1q+","Gain"),c("Xp-","Loss")));
capri.with.bootstrap.no.estimations = tronco.capri(data);
capri.without.bootstrap = tronco.capri(data,do.boot=FALSE,do.estimation=TRUE);

#perform the estimations given the error rates
estimated.error.rates.pf = list(error.fp=0,error.fn=0);
estimated.error.rates.bic = list(error.fp=0,error.fn=0);
capri.with.bootstrap.no.estimations.with.error.rates = tronco.estimation(capri.with.bootstrap.no.estimations,list(error.rates.pf=estimated.error.rates.pf,error.rates.bic=estimated.error.rates.bic));
#perform the estimations without giving the error rates
capri.with.bootstrap.no.estimations.without.error.rates = tronco.estimation(capri.with.bootstrap.no.estimations);

#perform the estimation by non-parametric and parametric bootstraps with CAPRI
set.seed("12345");
capri.without.bootstrap.non.parametric = tronco.bootstrap(capri.without.bootstrap);
set.seed("12345");
capri.without.bootstrap.parametric = tronco.bootstrap(capri.without.bootstrap,type="parametric");

#plot the results
tronco.plot(capri.without.bootstrap.non.parametric);
tronco.plot(capri.without.bootstrap.non.parametric,pf=TRUE);
