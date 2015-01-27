#file where the dataset is saved
file.dataset.ovarian = "data/ovarian-data.txt";

#set the working directory
work.dir = '~/Desktop/_tool_test';
setwd(work.dir);

#load the dataset and set all the values for colnames and rownames
genotypes = read.table(file.dataset.ovarian);
colnames(genotypes) = gsub("V", "gene ", colnames(genotypes));
rownames(genotypes) = paste("patient ",rownames(genotypes),sep="");
annotations = array("",c(ncol(genotypes),2));
colnames(annotations) = c("event","type");
rownames(annotations) = colnames(genotypes);
annotations[1,1] = "8q+";
annotations[2,1] = "3q+";
annotations[3,1] = "5q-";
annotations[4,1] = "4q-";
annotations[5,1] = "8p-";
annotations[6,1] = "1q+";
annotations[7,1] = "Xp-";
annotations[c(1,2,6),2] = "Gain";
annotations[c(3,4,5,7),2] = "Loss";
types = array("",c(length(unique(annotations[,"type"])),1));
colnames(types) = "color";
rownames(types) = unique(annotations[,"type"]);
types[1,1] = "Red";
types[2,1] = "Blue";

#create the input variable
data = list(genotypes = genotypes, annotations = annotations, types = types);

#load TRONCO package
invisible(sapply(list.files(pattern="[.]R$",path="R",full.names=TRUE),source));

#perform the reconstruction with CAPRESE with estimations
caprese = tronco.caprese(data,do.estimation=TRUE);

#perform the estimation by non-parametric and parametric bootstraps with CAPRESE
set.seed("12345");
caprese.non.parametric = tronco.bootstrap(caprese);
set.seed("12345");
caprese.parametric = tronco.bootstrap(caprese,type="parametric");

#perform the reconstruction with CAPRI using its default values
data = hypothesis.add(data,"H1",OR(XOR(c("8q+","Gain"),c("4q-","Loss")),AND(c("3q+","Gain"),c("5q-","Loss")),c("8p-","Loss"),c("1q+","Gain")),c("Xp-","Loss"));
data = hypothesis.add(data,"H2",AND(XOR(c("8q+","Gain"),c("4q-","Loss")),OR(c("3q+","Gain"),c("5q-","Loss")),c("8p-","Loss"),c("1q+","Gain")),c("Xp-","Loss"));
data = hypothesis.add(data,"H3",OR(XOR(c("8q+","Gain"),c("4q-","Loss")),c("8p-","Loss")),"*");
data = hypothesis.add(data,"H4",OR(AND(c("8q+","Gain"),c("4q-","Loss")),c("8p-","Loss")),c("1q+","Gain"),c("Xp-","Loss"));
capri.with.bootstrap = tronco.capri(data);
capri.without.bootstrap = tronco.capri(data,do.boot=FALSE,do.estimation=TRUE);

#perform the estimation by non-parametric and parametric bootstraps with CAPRI
set.seed("12345");
capri.without.bootstrap.non.parametric = tronco.bootstrap(capri.without.bootstrap);
set.seed("12345");
capri.without.bootstrap.parametric = tronco.bootstrap(capri.without.bootstrap,type="parametric");
