# File: bee_medip_overlap_bs.R
# Auth: u.niazi@imperial.ac.uk
# DESC: overlap of medip and bs-seq regions to test prediction accuracy
# Date: 26/11/2015


# header file
source('bee_data_header.R')

# select a scaffold for comparison to reduce memory requirements
csScaffold = 'NC_007070.3'

### data from macs
oGRmacs = f_oGRMACS2toGRanges(file.choose())

### rename sequences
s = gsub("gi\\|\\d+\\|ref\\|(.+)\\|", replacement = '\\1', x = as.character(seqnames(oGRmacs)), perl=T)
oGRmacs.s = GRanges(s, ranges(oGRmacs), mcols=mcols(oGRmacs))
oGRmacs.s = oGRmacs.s[seqnames(oGRmacs.s) == csScaffold]

##################################
#### medip data after strand assignment using dismiss
i = grep(csScaffold, x = as.character(seqnames(oGRmacs)), perl=T)
oGRdismiss = f_oGRSeparateStrands(oGRmacs[i], csFile = file.choose(), bPaired = T)

## rename sequences
s = gsub("gi\\|\\d+\\|ref\\|(.+)\\|", replacement = '\\1', x = as.character(seqnames(oGRdismiss)), perl=T)
oGRdismiss.s = GRanges(s, ranges(oGRdismiss), strand=strand(oGRdismiss), mcols=mcols(oGRdismiss))


##### repeat on bisulphite dataset
## choose the correct bismark methyl extractor files
## created previously using bee_dismiss_vs_bs.R
load(file='Bee_data/Objects/oGRbis.s')
# ## if we want to remove low count potentially noisy data
# ## model the distribution of the data to see how many times
# ## a particular C is methylated and if it appears to follow a poisson
# ## or negative binomial distribution, then use that to remove noisy bits
# # choose a cutoff by modelling the distribution shape
# w = countOverlaps(unique(oGRbis.s), oGRbis.s)
# w2 = log(w)
# r = range(w2)
# s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
# r[1] = floor(r[1])
# r[2] = ceiling(r[2])
# hist(w2, prob=T, breaks=s, main='distribution of 5mCs over nucleotide positions', 
#      xlab='log of counts', ylab='')
# r = round(r)
# dp = dpois(r[1]:r[2], lambda = mean(w2))
# dn = dnbinom(r[1]:r[2], size = mean(w2), mu = mean(w2))
# lines(r[1]:r[2], dp, col='red', type='b')
# lines(r[1]:r[2], dn, col='blue', type='b')
# legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
# c = qpois(0.05, mean(w2), lower.tail = F)
# points(c, 0, pch=20, col='red')
# ## we decide to use a poisson distribution
# oGRbis.s = f_oGRRemovePoissonNoiseFromBismark(oGRbis.s)
# gc()

## choose medip regions that overlap with bs-seq
f = overlapsAny(oGRdismiss.s, oGRbis.s, ignore.strand=T)
table(f)

oGRmedip = oGRdismiss.s[f]

## count the number of Cs on each side of the strands in each type of medip peak
## i.e. plus, minus and both (double) stranded peak classes
plus.ot = countOverlaps(oGRmedip[strand(oGRmedip) == '+'], oGRbis.s[strand(oGRbis.s) == '+'], ignore.strand=T)
plus.ob = countOverlaps(oGRmedip[strand(oGRmedip) == '+'], oGRbis.s[strand(oGRbis.s) == '-'], ignore.strand=T)

minus.ot = countOverlaps(oGRmedip[strand(oGRmedip) == '-'], oGRbis.s[strand(oGRbis.s) == '+'], ignore.strand=T)
minus.ob = countOverlaps(oGRmedip[strand(oGRmedip) == '-'], oGRbis.s[strand(oGRbis.s) == '-'], ignore.strand=T)

double.ot = countOverlaps(oGRmedip[strand(oGRmedip) == '*'], oGRbis.s[strand(oGRbis.s) == '+'], ignore.strand=T)
double.ob = countOverlaps(oGRmedip[strand(oGRmedip) == '*'], oGRbis.s[strand(oGRbis.s) == '-'], ignore.strand=T)


## continue from here with t.test paired



## calculate overlaps
### for each defined genomic feature, count how many peaks overlap
prom = overlapsAny(oGRprom, oGRbis.s)
ds = overlapsAny(oGRds, oGRbis.s)
genes = overlapsAny(oGRgene, oGRbis.s)
exon = overlapsAny(oGRexon, oGRbis.s)
intron = overlapsAny(oGRintron, oGRbis.s)

mDat = cbind(prom=as.numeric(table(prom)), exon=as.numeric(table(exon)), intron=as.numeric(table(intron)),
             ds=as.numeric(table(ds)))

rownames(mDat) = c('False', 'True')

# calculate confidence intervals/errors by simulating the data
ivProb.bis.s = mDat['True',] / rowSums(mDat)['True']

# simulate data as a sample from a multinomial distribution
# P(y | Theta)
mDat.bis.s = mDat['True',]
mSim.bis.s = t(rmultinom(n = 1000, size = 1000, prob = ivProb.bis.s))
# convert to probability scale for plotting
mSim.bis.s = mSim.bis.s / rowSums(mSim.bis.s)
# look at data distribution
l = barplot(colMeans(mSim.bis.s), ylim=c(0,0.7))
# calculate quantiles for drawing
m = apply(mSim.bis.s, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)

##### plot bs, medip and dismiss together
c = grey.colors(3)
mBar = rbind(ivProb.medip.s, ivProb.dismiss.s, ivProb.bis.s)
rownames(mBar) = c('MeDIP', 'DISMISS', 'BS-Seq')
l2 = barplot(mBar, beside=T, ylim=c(0, 0.7), col=c)
legend('topright', legend = rownames(mBar), fill = c)
l = l2[1,]
## make error bars
m = apply(mSim.medip.s, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)

l = l2[2,]
## make error bars
m = apply(mSim.dismiss.s, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)

l = l2[3,]
m = apply(mSim.bis.s, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)

############################################################
### count number of methylated 5mCs on each side of Dismiss peak
oGRbis.s = f_LoadObject(file.choose())
f = overlapsAny(oGRdismiss.s, oGRbis.s, ignore.strand=T)
oGRdismiss.s.bs = oGRdismiss.s[f]

# count how many 5mCs on from each strand
OT = countOverlaps(oGRdismiss.s.bs, (oGRbis.s[strand(oGRbis.s) == '+']), ignore.strand=T)
OB = countOverlaps(oGRdismiss.s.bs, (oGRbis.s[strand(oGRbis.s) == '-']), ignore.strand=T)
OT = log(OT+1)
OB = log(OB+1)
dfCounts = data.frame(OT, OB, strand=as.factor(strand(oGRdismiss.s.bs)))
dfCounts$Diff = dfCounts$OT - dfCounts$OB
## summarize the data
BS.Plus = tapply(dfCounts$OT, dfCounts$strand, mean)
BS.Minus = tapply(dfCounts$OB, dfCounts$strand, mean)
BS.Diff = tapply(dfCounts$Diff, dfCounts$strand, mean)
BS.Diff.SD = tapply(dfCounts$Diff, dfCounts$strand, mad)

## perform wilcoxon rank test and t test on data
wilcox.test(dfCounts$OT[dfCounts$strand=='*'], dfCounts$OB[dfCounts$strand=='*'], paired = T)$p.value
wilcox.test(dfCounts$OT[dfCounts$strand=='+'], dfCounts$OB[dfCounts$strand=='+'], paired = T)$p.value
wilcox.test(dfCounts$OT[dfCounts$strand=='-'], dfCounts$OB[dfCounts$strand=='-'], paired = T)$p.value

t.test(dfCounts$OT[dfCounts$strand=='*'], dfCounts$OB[dfCounts$strand=='*'], paired = T)$p.value
t.test(dfCounts$OT[dfCounts$strand=='+'], dfCounts$OB[dfCounts$strand=='+'], paired = T)$p.value
t.test(dfCounts$OT[dfCounts$strand=='-'], dfCounts$OB[dfCounts$strand=='-'], paired = T)$p.value

x.plus = rnorm(1000, BS.Diff['+'], BS.Diff.SD['+'])
x.minus = rnorm(1000, BS.Diff['-'], BS.Diff.SD['-'])
x.ds = rnorm(1000, BS.Diff['*'], BS.Diff.SD['*'])
mSim = cbind(x.plus, x.minus, x.ds)
### plot the data
names(BS.Diff) = c('Plus', 'Minus', 'Double')
l = barplot(BS.Diff, ylim=c(-1, 1))

## make error bars
m = apply(mSim, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)


