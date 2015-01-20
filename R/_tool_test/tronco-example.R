#file where the dataset is saved
file.dataset.ovarian = "data/ovarian-data.txt";

#set the working directory
work.dir = '~/Desktop/_tool_test';
setwd(work.dir);

#load the dataset and set all the values for colnames and rownames
genotypes = read.table(file.dataset.ovarian);
colnames(genotypes) = gsub("V", "Gene ", colnames(genotypes));
rownames(genotypes) = paste("Patient ",rownames(genotypes),sep="");
annotations = array("",c(ncol(genotypes),2));
colnames(annotations) = c("Event","Type");
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
types = array("",c(length(unique(annotations[,"Type"])),1));
colnames(types) = "Color";
rownames(types) = unique(annotations[,"Type"]);
types[1,1] = "Red";
types[2,1] = "Blue";

#create the input variable
data = list(genotypes = genotypes, annotations = annotations, types = types);

#load TRONCO package
invisible(sapply(list.files(pattern="[.]R$",path="R",full.names=TRUE),source));

#perform the reconstruction with CAPRESE
caprese = caprese.fit(data$genotypes);

#perform the reconstruction with CAPRI
my.hypotheses = hypothesis.add(data$genotypes,"H1",OR(XOR("Gene 1","Gene 4"),AND("Gene 2","Gene 3"),"Gene 5","Gene 6"),"Gene 7");
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H2",AND(XOR("Gene 1","Gene 4"),OR("Gene 2","Gene 3"),"Gene 5","Gene 6"),"Gene 7",my.hypotheses$hypotheses);
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H3",OR(XOR("Gene 1","Gene 4"),"Gene 5"),"*",my.hypotheses$hypotheses);
capri = capri.fit(my.hypotheses$dataset,my.hypotheses$hypotheses);
