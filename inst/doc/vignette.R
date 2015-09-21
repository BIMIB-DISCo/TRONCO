### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: vignette.Rnw:85-88
###################################################
library(TRONCO)
data(aCML)
hide.progress.bar <<- TRUE


###################################################
### code chunk number 3: vignette.Rnw:93-94
###################################################
show(aCML)


###################################################
### code chunk number 4: vignette.Rnw:99-100
###################################################
as.events(aCML)


###################################################
### code chunk number 5: vignette.Rnw:105-106
###################################################
as.genes(aCML)


###################################################
### code chunk number 6: vignette.Rnw:111-112
###################################################
as.gene(aCML, genes='SETBP1')


###################################################
### code chunk number 7: vignette.Rnw:117-118
###################################################
gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')


###################################################
### code chunk number 8: vignette.Rnw:123-124
###################################################
alterations = events.selection(as.alterations(aCML), filter.freq = .05)


###################################################
### code chunk number 9: vignette.Rnw:129-130
###################################################
dummy = oncoprint(alterations,font.row=12,cellheight=20,cellwidth=4)


###################################################
### code chunk number 10: vignette.Rnw:132-133
###################################################
capture.output(oncoprint(alterations,font.row=12,cellheight=20,cellwidth=4, file='onco-1.pdf'), file='NUL')


###################################################
### code chunk number 11: vignette.Rnw:143-145
###################################################
hypo = events.selection(aCML, filter.in.names=c(as.genes(alterations), gene.hypotheses))
hypo = annotate.description(hypo, 'CAPRI - Bionformatics aCML data (selected events)')


###################################################
### code chunk number 12: vignette.Rnw:150-152
###################################################
dummy = oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, 
	font.row=12, font.column=5, cellheight=20, cellwidth=4)


###################################################
### code chunk number 13: vignette.Rnw:154-155
###################################################
capture.output(oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, font.row=12, font.column=5, cellheight=20, cellwidth=4, file='onco-2.pdf'), file='NUL')


###################################################
### code chunk number 14: vignette.Rnw:163-164
###################################################
hypo = hypothesis.add(hypo, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))


###################################################
### code chunk number 15: vignette.Rnw:169-170 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'NRAS or KRAS',  OR('NRAS', 'KRAS'))


###################################################
### code chunk number 16: vignette.Rnw:175-177
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')),
	font.row=12,cellheight=20, cellwidth=4)


###################################################
### code chunk number 17: vignette.Rnw:179-180
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')),font.row=12,cellheight=20, cellwidth=4, file='onco-3.pdf'), file='NUL')


###################################################
### code chunk number 18: vignette.Rnw:188-189
###################################################
hypo = hypothesis.add(hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 19: vignette.Rnw:191-192 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'SF3B1 or ASXL1', OR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 20: vignette.Rnw:197-200
###################################################
as.events(hypo, genes = 'TET2') 
hypo = hypothesis.add(hypo, 'TET2 xor IDH2', XOR('TET2', 'IDH2'), '*')
hypo = hypothesis.add(hypo,  'TET2 or IDH2', OR('TET2', 'IDH2'), '*')


###################################################
### code chunk number 21: vignette.Rnw:203-205
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')),
	font.row=12, cellheight=20,cellwidth=4)


###################################################
### code chunk number 22: vignette.Rnw:207-208
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')),font.row=12, cellheight=20,cellwidth=4, file='onco-4.pdf'), file='NUL')


###################################################
### code chunk number 23: vignette.Rnw:216-217
###################################################
hypo = hypothesis.add.homologous(hypo)


###################################################
### code chunk number 24: vignette.Rnw:222-224
###################################################
dummy = oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, 
	font.row=10, font.column=5, cellheight=15, cellwidth=4)


###################################################
### code chunk number 25: vignette.Rnw:226-227
###################################################
capture.output(oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, font.row=10, font.column=5, cellheight=15, cellwidth=4, file='onco-5.pdf'), file='NUL')


###################################################
### code chunk number 26: vignette.Rnw:237-238
###################################################
model = tronco.capri(hypo, boot.seed = 12345, nboot=10)


###################################################
### code chunk number 27: figplot
###################################################
tronco.plot(model, fontsize = 13, scale.nodes = .6, regularization="bic", 
	confidence = c('tp', 'pr', 'hg'), height.logic = 0.25, legend.cex = .5, 
	pathways =  list(priors= gene.hypotheses), label.edge.size=5)



###################################################
### code chunk number 28: vignette.Rnw:258-259
###################################################
model.boot = tronco.bootstrap(model, nboot=10)


###################################################
### code chunk number 29: figplotboot
###################################################
tronco.plot(model.boot, fontsize = 13, scale.nodes = .6, regularization="bic", 
	confidence=c('npb'), height.logic = 0.25, legend.cex = .5, 
	pathways = list(priors= gene.hypotheses), label.edge.size=10)



