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
my.hypotheses = hypothesis.add(dataset,"H1",OR(XOR("V1","V4"),AND("V2","V3"),"V5","V6"),"V7");
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H2",AND(XOR("V1","V4"),OR("V2","V3"),"V5","V6"),"V7",my.hypotheses$hypotheses);
my.hypotheses = hypothesis.add(my.hypotheses$dataset,"H3",OR(XOR("V1","V4"),"V5"),"*",my.hypotheses$hypotheses);
capri = capri.fit(my.hypotheses$dataset,my.hypotheses$hypotheses);

#print the lifted adjacency matrix for each hypothesis
print(my.hypotheses$hypotheses$hstructure[["H1"]]);
print(my.hypotheses$hypotheses$hstructure[["H2"]]);
print(my.hypotheses$hypotheses$hstructure[["H3"]]);

# set colnames e rownames to adjacency matrix
colnames(capri$adj.matrix$adj.matrix.bic) = colnames(my.hypotheses$data)
rownames(capri$adj.matrix$adj.matrix.bic) = colnames(my.hypotheses$data)

# load hypotheses.
source('../hypotheses.expansion.R')
source('../tronco.plot.R')


hypo_mat = hypotheses.expansion(capri$adj.matrix$adj.matrix.bic, 
                                  my.hypotheses$hypotheses$num.hypotheses, 
                                  my.hypotheses$hypotheses$hstructure)

hypo_graph = graph.adjacency(hypo_mat)
graph <- igraph.to.graphNEL(hypo_graph)
graph <- layoutGraph(graph)
graph.par(list(nodes=list(fontsize=70, textCol="black")))
renderGraph(graph)
