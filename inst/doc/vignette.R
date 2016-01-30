### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: vignette.Rnw:82-85
###################################################
library(TRONCO)
data(aCML)
hide.progress.bar <<- TRUE


###################################################
### code chunk number 3: vignette.Rnw:90-91
###################################################
view(aCML)


###################################################
### code chunk number 4: vignette.Rnw:96-97
###################################################
as.events(aCML)


###################################################
### code chunk number 5: vignette.Rnw:102-103
###################################################
as.genes(aCML)


###################################################
### code chunk number 6: vignette.Rnw:108-109
###################################################
as.gene(aCML, genes='SETBP1')


###################################################
### code chunk number 7: vignette.Rnw:114-115
###################################################
gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')


###################################################
### code chunk number 8: vignette.Rnw:120-121
###################################################
alterations = events.selection(as.alterations(aCML), filter.freq = .05)


###################################################
### code chunk number 9: vignette.Rnw:126-127
###################################################
dummy = oncoprint(alterations,font.row=12,cellheight=20,cellwidth=4)


###################################################
### code chunk number 10: vignette.Rnw:129-130
###################################################
capture.output(oncoprint(alterations,font.row=12,cellheight=20,cellwidth=4, file='onco-1.pdf'), file='NUL')


###################################################
### code chunk number 11: vignette.Rnw:140-142
###################################################
hypo = events.selection(aCML, filter.in.names=c(as.genes(alterations), gene.hypotheses))
hypo = annotate.description(hypo, 'CAPRI - Bionformatics aCML data (selected events)')


###################################################
### code chunk number 12: vignette.Rnw:147-149
###################################################
dummy = oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, 
	font.row=12, font.column=5, cellheight=20, cellwidth=4)


###################################################
### code chunk number 13: vignette.Rnw:151-152
###################################################
capture.output(oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, font.row=12, font.column=5, cellheight=20, cellwidth=4, file='onco-2.pdf'), file='NUL')


###################################################
### code chunk number 14: vignette.Rnw:160-161
###################################################
hypo = hypothesis.add(hypo, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))


###################################################
### code chunk number 15: vignette.Rnw:166-167 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'NRAS or KRAS',  OR('NRAS', 'KRAS'))


###################################################
### code chunk number 16: vignette.Rnw:172-174
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')),
	font.row=12,cellheight=20, cellwidth=4)


###################################################
### code chunk number 17: vignette.Rnw:176-177
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')),font.row=12,cellheight=20, cellwidth=4, file='onco-3.pdf'), file='NUL')


###################################################
### code chunk number 18: vignette.Rnw:185-186
###################################################
hypo = hypothesis.add(hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 19: vignette.Rnw:188-189 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'SF3B1 or ASXL1', OR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 20: vignette.Rnw:194-197
###################################################
as.events(hypo, genes = 'TET2') 
hypo = hypothesis.add(hypo, 'TET2 xor IDH2', XOR('TET2', 'IDH2'), '*')
hypo = hypothesis.add(hypo,  'TET2 or IDH2', OR('TET2', 'IDH2'), '*')


###################################################
### code chunk number 21: vignette.Rnw:200-202
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')),
	font.row=12, cellheight=20,cellwidth=4)


###################################################
### code chunk number 22: vignette.Rnw:204-205
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')),font.row=12, cellheight=20,cellwidth=4, file='onco-4.pdf'), file='NUL')


###################################################
### code chunk number 23: vignette.Rnw:213-214
###################################################
hypo = hypothesis.add.homologous(hypo)


###################################################
### code chunk number 24: vignette.Rnw:219-221
###################################################
dummy = oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, 
	font.row=10, font.column=5, cellheight=15, cellwidth=4)


###################################################
### code chunk number 25: vignette.Rnw:223-224
###################################################
capture.output(oncoprint(hypo, gene.annot = list(priors= gene.hypotheses), sample.id = T, font.row=10, font.column=5, cellheight=15, cellwidth=4, file='onco-5.pdf'), file='NUL')


###################################################
### code chunk number 26: vignette.Rnw:234-235
###################################################
model = tronco.capri(hypo, boot.seed = 12345, nboot=10)


###################################################
### code chunk number 27: figplot
###################################################
tronco.plot(model, fontsize = 13, scale.nodes = .6, regularization="bic", 
	confidence = c('tp', 'pr', 'hg'), height.logic = 0.25, legend.cex = .5, 
	pathways =  list(priors= gene.hypotheses), label.edge.size=5)



###################################################
### code chunk number 28: vignette.Rnw:255-256
###################################################
model.boot = tronco.bootstrap(model, nboot=10)


###################################################
### code chunk number 29: figplotboot
###################################################
tronco.plot(model.boot, fontsize = 13, scale.nodes = .6, regularization="bic", 
	confidence=c('npb'), height.logic = 0.25, legend.cex = .5, 
	pathways = list(priors= gene.hypotheses), label.edge.size=10)



