### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: vignette.Rnw:96-99
###################################################
library(TRONCO)
data(aCML)
hide.progress.bar <<- TRUE


###################################################
### code chunk number 3: vignette.Rnw:103-104
###################################################
show(data)


###################################################
### code chunk number 4: vignette.Rnw:108-109
###################################################
as.events(data)


###################################################
### code chunk number 5: vignette.Rnw:113-114
###################################################
as.genes(data)


###################################################
### code chunk number 6: vignette.Rnw:118-119
###################################################
as.gene(data, genes='SETBP1')


###################################################
### code chunk number 7: vignette.Rnw:124-125
###################################################
gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')


###################################################
### code chunk number 8: vignette.Rnw:132-133
###################################################
alterations = events.selection(as.alterations(data), filter.freq = .05)


###################################################
### code chunk number 9: onco1
###################################################
dummy = oncoprint(alterations)


###################################################
### code chunk number 10: vignette.Rnw:151-153
###################################################
hypo = events.selection(data, filter.in.names=c(as.genes(alterations), gene.hypotheses))
hypo = annotate.description(hypo, 'CAPRI - Bionformatics aCML data (selected events)')


###################################################
### code chunk number 11: onco2
###################################################
dummy = oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T)


###################################################
### code chunk number 12: vignette.Rnw:169-170
###################################################
hypo = hypothesis.add(hypo, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))


###################################################
### code chunk number 13: vignette.Rnw:175-178
###################################################
### do not run
# hypo = hypothesis.add(hypo, 'NRAS or KRAS',  OR('NRAS', 'KRAS'))
###


###################################################
### code chunk number 14: onco3
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('KRAS', 'NRAS')))


###################################################
### code chunk number 15: vignette.Rnw:191-195
###################################################
hypo = hypothesis.add(hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', OR('ASXL1')), '*')
### do not run
# hypo = hypothesis.add(hypo, 'SF3B1 or ASXL1', OR('SF3B1', OR('ASXL1')), '*')
###


###################################################
### code chunk number 16: vignette.Rnw:202-207
###################################################
as.events(hypo, genes = 'TET2') 
hypo = hypothesis.add(hypo, 'TET2 xor IDH2', XOR('TET2', 'IDH2'), '*')
### do not run
# hypo = hypothesis.add(hypo,  'TET2 or IDH2', OR('TET2', 'IDH2'), '*'))
###


###################################################
### code chunk number 17: onco4
###################################################
dummy = oncoprint(events.selection(hypo, filter.in.names = c('TET2', 'IDH2')))


###################################################
### code chunk number 18: vignette.Rnw:219-220
###################################################
hypo = hypothesis.add.homologous(hypo)


###################################################
### code chunk number 19: onco5
###################################################
dummy = oncoprint(hypo,  gene.annot = list(priors= gene.hypotheses), sample.id = T)


###################################################
### code chunk number 20: vignette.Rnw:235-236
###################################################
model = tronco.capri(hypo, boot.seed = 12345, regularization='bic')


###################################################
### code chunk number 21: figplot
###################################################
tronco.plot(model,  
  fontsize = 13,  
  scale.nodes = .6,   
  confidence = c('tp', 'pr', 'hg'),   
  height.logic = 0.25,  
  legend.cex = .5, 
  pathways =  list(priors= gene.hypotheses))


###################################################
### code chunk number 22: vignette.Rnw:258-259
###################################################
model.boot = tronco.bootstrap(model, nboot=6)


###################################################
### code chunk number 23: figplotboot
###################################################
tronco.plot(model.boot,
            fontsize = 13,  
            scale.nodes = .6,
            confidence=c('npb'),
            height.logic = 0.25,  
            legend.cex = .5)


