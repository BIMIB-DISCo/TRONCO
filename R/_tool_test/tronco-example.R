#dataset
file.dataset.ovarian = "data/ovarian-data.txt";

#set the working directory
work.dir = '~/Desktop/_tool_test';
setwd(work.dir);

#load the dataset
dataset = read.table(file.dataset.ovarian);

#load TRONCO package
invisible(sapply(list.files(pattern="[.]R$",path="R",full.names=TRUE),source));

#perform the reconstruction with CAPRESE
caprese = caprese.fit(dataset);

#perform the reconstruction with CAPRI
my.hypotheses = hypothesis.add(dataset,"H1",OR(dataset,XOR(dataset,"V1","V4"),AND(dataset,"V2","V3"),"V5","V6"),"V7");
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H2",AND(my.hypotheses$dataset,XOR(my.hypotheses$dataset,"V1","V4"),OR(my.hypotheses$dataset,"V2","V3"),"V5","V6"),"V7",hypotheses=my.hypotheses$hypotheses);
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H3",OR(my.hypotheses$dataset,XOR(my.hypotheses$dataset,"V1","V4"),"V5"),"*",hypotheses=my.hypotheses$hypotheses);
capri = capri.fit(my.hypotheses$dataset,my.hypotheses$hypotheses);
