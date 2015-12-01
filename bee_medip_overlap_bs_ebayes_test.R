# File: bee_medip_overlap_bs.R
# Auth: u.niazi@imperial.ac.uk
# DESC: overlap of medip and bs-seq regions to test prediction accuracy
# Date: 26/11/2015


# header file
source('bee_data_header.R')

# select a scaffold for comparison to reduce memory requirements
csScaffold = 'NC_007070.3'

## choose the appropriate title and data
csPlot.title = 'SRR850130 - MeDIP Apis millefera ligustica nurse honeybee'
csPlot.title = 'SRR850131 - MeDIP Apis millefera ligustica forager honeybee'
csPlot.title = 'SRR850132 - MeDIP Apis millefera ligustica reverted nurse honeybee'

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

## if we want to remove low count potentially noisy data
## model the distribution of the data to see how many times
## a particular C is methylated and if it appears to follow a poisson
## or negative binomial distribution, then use that to remove noisy bits
# choose a cutoff by modelling the distribution shape
w = countOverlaps(unique(oGRbis.s), oGRbis.s)
w2 = log(w)
r = range(w2)
s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
r[1] = floor(r[1])
r[2] = ceiling(r[2])
hist(w2, prob=T, breaks=s, main='distribution of 5mCs over nucleotide positions', 
     xlab='log of counts', ylab='')
r = round(r)
dp = dpois(r[1]:r[2], lambda = mean(w2))
dn = dnbinom(r[1]:r[2], size = mean(w2), mu = mean(w2))
lines(r[1]:r[2], dp, col='red', type='b')
lines(r[1]:r[2], dn, col='blue', type='b')
legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
c = qpois(0.05, mean(w2), lower.tail = F)
points(c, 0, pch=20, col='red')
## we decide to use a poisson distribution
oGRbis.s = f_oGRRemovePoissonNoiseFromBismark(oGRbis.s)
gc()

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

signif(t.test(plus.ot, plus.ob, paired=T)$p.value, 3)
signif(t.test(minus.ot, minus.ob, paired=T)$p.value, 3)
signif(t.test(double.ot, double.ob, paired=T)$p.value, 3)

## bar plots of proprtions
getPost = function(s, f){
  # get the estimate for the prior using e-bayes to estimate hyperparameters
  t = s + f
  s.p = s/t
  s.p.mean = mean(s.p)
  s.p.var = var(s.p)
  al.be = getalphabeta(s.p.mean, s.p.var)
  # simulate the theta for each data using the hyperparameters
#   prior = rbeta(10000, al.be$alpha, al.be$beta)
#   y = round(mean(s), 0)
#   n = round(mean(s) + mean(f), 0)
#   lik = sapply(prior, function(x) dbinom(y, n, x))
#   rs = sample(prior, size = 1000, replace = T, prob = lik)
  rs = rbeta(1000, mean(s)+al.be$alpha, mean(f)+al.be$beta)
  rf = 1-rs #rbeta(1000, f+2, s+2)
  return(list(s=rs, f=rf))  
}

getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(list(alpha=al, beta=be))
}

plus = getPost((plus.ot+1), (plus.ob+1))
minus = getPost((minus.ot+1), (minus.ob+1))
double = getPost((double.ot+1), (double.ob+1))

Plus.Stranded = sapply(plus, median)
Minus.Stranded = sapply(minus, median)
Double.Stranded = sapply(double, median)

mDat = cbind(Plus.Stranded, Minus.Stranded, Double.Stranded)
rownames(mDat) = c('mCPlus', 'mCMinus')

col = grey.colors(2)
l = barplot(mDat, beside=T, col=col, main=csPlot.title, ylim=c(0,1),
            horiz =F, cex.names=0.7)

f_barplot_errorbars = function(x.loc, y.loc, ...){
  segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
  segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
  segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
}

sapply(seq_along(plus), function(x) f_barplot_errorbars(l[x,1], quantile(plus[[x]], c(0.025, 0.975))))
sapply(seq_along(minus), function(x) f_barplot_errorbars(l[x,2], quantile(minus[[x]], c(0.025, 0.975))))
sapply(seq_along(double), function(x) f_barplot_errorbars(l[x,3], quantile(double[[x]], c(0.025, 0.975))))



