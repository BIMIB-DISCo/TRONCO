
source('../selection/events.selection.R')
source('../visualization/oncoprint.R')

# Load data - makes 2 variables available
# tcga: whole dataset
# tcga.sub1: subtype 1, whole dataset
source('./example-dataloading.R')

# Events selection: examples
out.names = c('TTN', 'ADAM6')
in.names = c('PTEN', 'TP53', 'KRAS', 'ARID1A')
minimum.freq = .2 # minimum mutation frequency

tcga.sel = events.selection(tcga, minimum.freq, in.names, out.names)
oncoprint(tcga.sel)

tcga.sub1.sel = events.selection(tcga.sub1, minimum.freq, in.names, out.names)
oncoprint(tcga.sub1.sel)

tcga.sub2.sel = events.selection(tcga.sub2, minimum.freq, in.names, out.names)
oncoprint(tcga.sub2.sel)
