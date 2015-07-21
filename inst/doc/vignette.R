### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: vignette.Rnw:99-102
###################################################
library(TRONCO)
data(aCML)
hide.progress.bar <<- TRUE


###################################################
### code chunk number 3: vignette.Rnw:106-107
###################################################
show(aCML)


###################################################
### code chunk number 4: vignette.Rnw:111-112
###################################################
as.events(aCML)


###################################################
### code chunk number 5: vignette.Rnw:116-117
###################################################
as.genes(aCML)


###################################################
### code chunk number 6: vignette.Rnw:121-122
###################################################
as.gene(aCML, genes='SETBP1')


###################################################
### code chunk number 7: vignette.Rnw:127-128
###################################################
gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')


###################################################
### code chunk number 8: vignette.Rnw:135-136
###################################################
alterations = events.selection(as.alterations(aCML), filter.freq = .05)


###################################################
### code chunk number 9: vignette.Rnw:141-142
###################################################
dummy = oncoprint(alterations)


###################################################
### code chunk number 10: vignette.Rnw:144-145
###################################################
capture.output(oncoprint(alterations, file='onco-1.pdf'), file='NUL')


###################################################
### code chunk number 11: vignette.Rnw:157-159
###################################################
hypo = events.selection(aCML, filter.in.names=c(as.genes(alterations), gene.hypotheses))
hypo = annotate.description(hypo, 'CAPRI - Bionformatics aCML data (selected events)')


###################################################
### code chunk number 12: vignette.Rnw:164-165
###################################################
dummy = oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T)


###################################################
### code chunk number 13: vignette.Rnw:167-168
###################################################
capture.output(oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T, file='onco-2.pdf'), file='NUL')


###################################################
### code chunk number 14: vignette.Rnw:178-179
###################################################
hypo = hypothesis.add(hypo, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))


###################################################
### code chunk number 15: vignette.Rnw:184-185 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'NRAS or KRAS',  OR('NRAS', 'KRAS'))


###################################################
### code chunk number 16: vignette.Rnw:190-191
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')))


###################################################
### code chunk number 17: vignette.Rnw:193-194
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')), file='onco-3.pdf'), file='NUL')


###################################################
### code chunk number 18: vignette.Rnw:201-202
###################################################
hypo = hypothesis.add(hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 19: vignette.Rnw:204-205 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo, 'SF3B1 or ASXL1', OR('SF3B1', OR('ASXL1')), '*')


###################################################
### code chunk number 20: vignette.Rnw:212-214
###################################################
as.events(hypo, genes = 'TET2') 
hypo = hypothesis.add(hypo, 'TET2 xor IDH2', XOR('TET2', 'IDH2'), '*')


###################################################
### code chunk number 21: vignette.Rnw:216-217 (eval = FALSE)
###################################################
## hypo = hypothesis.add(hypo,  'TET2 or IDH2', OR('TET2', 'IDH2'), '*')


###################################################
### code chunk number 22: vignette.Rnw:220-221
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')))


###################################################
### code chunk number 23: vignette.Rnw:223-224
###################################################
capture.output(oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')), file='onco-4.pdf'), file='NUL')


###################################################
### code chunk number 24: vignette.Rnw:232-233
###################################################
hypo = hypothesis.add.homologous(hypo)


###################################################
### code chunk number 25: vignette.Rnw:237-238
###################################################
dummy = oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T)


###################################################
### code chunk number 26: vignette.Rnw:240-241
###################################################
capture.output(oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T, file='onco-5.pdf'), file='NUL')


###################################################
### code chunk number 27: vignette.Rnw:251-252
###################################################
model = tronco.capri(hypo, boot.seed = 12345, regularization='bic', nboot=6)


###################################################
### code chunk number 28: figplot
###################################################
tronco.plot(model,  
  fontsize = 13,  
  scale.nodes = .6,   
  confidence = c('tp', 'pr', 'hg'),   
  height.logic = 0.25,  
  legend.cex = .5, 
  pathways =  list(priors= gene.hypotheses))


###################################################
### code chunk number 29: vignette.Rnw:274-275
###################################################
model.boot = tronco.bootstrap(model, nboot=6)


###################################################
### code chunk number 30: figplotboot
###################################################
tronco.plot(model.boot,
            fontsize = 13,  
            scale.nodes = .6,
            confidence=c('npb'),
            height.logic = 0.25,  
            legend.cex = .5)


