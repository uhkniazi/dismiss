# File: reviewers.R
# Auth: u.niazi@imperial.ac.uk
# DESC: analyses to answer questions from reviewers
# Date: 26/04/2016


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


############## Reviewer 1 Questions:
# An illustrative example of the distribution of plus/minus reads and the computation 
# of the different λs from that would be helpful
Y.p = 60
Y.m = 40
lambda = 40:60
lik = function(lam, dat) dpois(dat, lam)
lik.p = lik(lambda, Y.p)
lik.m = lik(lambda, Y.m)
lik.p = lik.p/max(lik.p)
lik.m = lik.m/max(lik.m)
plot(lambda, lik.p, type='l', col=1, lwd=2, xlab=expression(lambda), 
     ylab=expression(alpha), main='Maximum Likelihood Estimate',
     sub='Plausible Values of Lambda')
lines(lambda, lik.m, col=2, lwd=2)
legend('topleft', legend = c('Plus', 'Minus'), fill=c(1, 2))
abline(h = 0.1, lty=2)
# maximum liklihood for each data
lambda[which.max(lik.p)]
lambda[which.max(lik.m)]
sp = lambda[which(lik.p >= 0.1)]
sm = lambda[which(lik.m >= 0.1)]
table(sp %in% sm)

# What is the average number of reads (plus and minus) per peak 
# in MeDIP-Seq analysis and how does this number influence the algorithm?
# Can you give a histogram for all peak regions.
ivPlus = oGRdismiss.s$mcols.mp
ivMinus = oGRdismiss.s$mcols.mm
# trim the last quantile
fp = ivPlus < quantile(ivPlus, 0.95)
fm = ivMinus < quantile(ivMinus, 0.95)
# merge the boolean variable
fb = fp & fm
# plots
par(mfrow=c(1,2))
hist(ivPlus[fb], xlab='Read Depth', main='Plus Strand Read Depth', prob=T)
hist(ivMinus[fb], xlab='Read Depth', main='Minus Strand Read Depth', prob=T)
par(mfrow=c(1,2))
# plot(ivPlus[fb], ivMinus[fb], pch=20, cex=0.5, col=oGRdismiss.s$mcols.strand_fac[fb]+1,
#      xlab='Plus Stranded Reads', ylab='Minus Stranded reads')
# #     main='Number of Reads from Plus and Minus Strands at each Peak Position')
# plot.new()
# legend('center', legend = c('Plus', 'Minus', 'Double'), fill=c('red', 'black', 'green'))
grey.col = grey.colors(3, start = 0.1, end=0.7)
grey.col = grey.col[c(2,3,1)]
col.v = grey.col[oGRdismiss.s$mcols.strand_fac[fb]+1]
plot(ivPlus[fb], ivMinus[fb], pch=20, cex=0.5, col=col.v,
     xlab='Plus Stranded Reads', ylab='Minus Stranded reads')
#     main='Number of Reads from Plus and Minus Strands at each Peak Position')
plot.new()
legend('center', legend = c('Minus', 'Plus', 'Double'), fill=grey.col)
summary(ivPlus[fb])
summary(ivMinus[fb])
summary(ivPlus[fb]+ivMinus[fb])
# cut read depth into quantiles
rd = ivPlus[fb] + ivMinus[fb]
grps = cut(rd, breaks = quantile(rd, 0:20/20), include.lowest = T)
fStrands = strand(oGRdismiss.s[fb])
# strands in each quantile
lStrands = tapply(fStrands, grps, table)
mStrands = do.call(cbind, lStrands)
mStrands = mStrands/sum(rowSums(mStrands))
# P(Quantile)
pq = colSums(mStrands)
# P(strand|Quantile)
mStrands.cond = sweep(mStrands, 2, pq, FUN = '/')
matplot(t(mStrands.cond[,-20]), type = 'l', lty = 1, lwd = 2, col=1:3, 
        xaxt='n', xlab='Quantiles of Read Depth', ylab='Proportion of Strand Distribution',
        main='Relationship between Strand Distribution and Read Depth')
axis(1, at = 1:20, labels = names(quantile(rd, 0:20/20))[2:21], las=2)


# Peak regions are between 260 and 2500bp. It is likely that larger regions 
# contain both CpG and non-CpG sites and thus, both symmetric and asymmetric
# methylation might occur and with MeDIP this cannot be distinguished. 
# How is this reflected in the different statistics. Could that add some bias?
rd = width(oGRdismiss.s)
grps = cut(rd, breaks = quantile(rd, probs = c(0, 0.5, 0.95, 1)), include.lowest = T)
fStrands = strand(oGRdismiss.s)
# strands in each quantile
lStrands = tapply(fStrands, grps, table)
mStrands = do.call(cbind, lStrands)
chisq.test(mStrands)
st = colSums(mStrands[1:2,])
ns = mStrands[3,]
mStrands = rbind(st, ns)
rownames(mStrands) = c('Stranded', 'Non-Stranded')
mStrands = mStrands/sum(rowSums(mStrands))
# P(Quantile)
pq = colSums(mStrands)
# P(strand|Quantile)
mStrands.cond = sweep(mStrands, 2, pq, FUN = '/')
matplot(t(mStrands.cond), type = 'l', lty = 1, lwd = 2, col=1:2, 
        xaxt='n', xlab='Quantiles of Peak Width', ylab='Proportion of Strand Distribution',
        main='Relationship between Strand Distribution and Peak Width')
axis(1, at = 1:3, labels = c('[0-50]', '(50-95]', '(95-100]'))







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



