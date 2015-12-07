# File: bee_medip_overlap_bs.R
# Auth: u.niazi@imperial.ac.uk
# DESC: overlap of medip and bs-seq regions to test prediction accuracy
# Date: 26/11/2015


# header file
source('bee_data_header.R')

# select a scaffold for comparison to reduce memory requirements
csScaffold = 'NC_007070.3'
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
  # simulate the theta for each data using a uniform beta prior
  rs = rbeta(1000, s+2, f+2)
  rf = rbeta(1000, f+2, s+2)
  return(list(s=rs, f=rf))  
}


plus = getPost(mean(plus.ot), mean(plus.ob))
minus = getPost(mean(minus.ot), mean(minus.ob))
double = getPost(mean(double.ot), mean(double.ob))

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



