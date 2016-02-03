# File: bee_dismiss_medip_bs_rates.R
# Auth: u.niazi@imperial.ac.uk
# DESC: data comparison of medip and bs-seq
# Date: 02/02/2016


# header file
source('bee_data_header.R')
jef.prior = c(alpha=0.5, beta=0.0001)
bs.cut = 3

# load gff file
csFile = 'Bee_data/Bee_annotations/Amel_4.5_chr_bee/ref_Amel_4.5_top_level.gff3'
oGRgff = import.gff3(csFile)
# select a scaffold for comparison to reduce memory requirements
csScaffold = 'NC_007070.3'
oGRgff = oGRgff[seqnames(oGRgff) %in% csScaffold]
oGRgene = oGRgff[oGRgff$type=='gene']
oGRexon = oGRgff[oGRgff$type=='exon']

### data from macs2
oGRmacs = f_oGRMACS2toGRanges(file.choose())

### rename sequences
s = gsub("gi\\|\\d+\\|ref\\|(.+)\\|", replacement = '\\1', x = as.character(seqnames(oGRmacs)), perl=T)
oGRmacs.s = GRanges(s, ranges(oGRmacs), mcols=mcols(oGRmacs))
oGRmacs.s = oGRmacs.s[seqnames(oGRmacs.s) == csScaffold]

### gene scaffold for overlaps
oGRgene = reduce(oGRgene)
oGRexon = reduce(oGRexon)
oGRintron = setdiff(oGRgene, oGRexon)
oGRintron = reduce(oGRintron)
oGRprom = flank(oGRgene, width = 2000, start = T)
oGRds = flank(oGRgene, width = 2000, start = F)
# any promoters or downstream regoins at ends of chromosomes
f = width(oGRprom) < 2000
oGRprom = oGRprom[!f]
f = width(oGRds) < 2000
oGRds = oGRds[!f]

# remove promoters that run into a gene previous gene
f = overlapsAny(oGRprom, oGRgene)
oGRprom = oGRprom[!f]
# remove downstream regions running into the next gene
f = overlapsAny(oGRds, oGRgene)
oGRds = oGRds[!f]

### for each defined genomic feature, count how many peaks overlap
### with medip (macs) data
prom = overlapsAny(oGRprom, oGRmacs.s)
ds = overlapsAny(oGRds, oGRmacs.s)
genes = overlapsAny(oGRgene, oGRmacs.s)
exon = overlapsAny(oGRexon, oGRmacs.s)
intron = overlapsAny(oGRintron, oGRmacs.s)

mDat = cbind(prom=as.numeric(table(prom)), exon=as.numeric(table(exon)), intron=as.numeric(table(intron)),
             ds=as.numeric(table(ds)))

rownames(mDat) = c('False', 'True')

ivDat = colSums(mDat)
ivProb.medip.s = sweep(mDat, 2, ivDat, '/')['True',] * 1000
# calculate posterior for theta i.e. rate using jeffereys non informative prior
post = sapply(ivProb.medip.s, function(x) {
  r = rgamma(1000, shape = jef.prior['alpha'] + sum(x) , jef.prior['beta']+ length(x))
  return(r)
})

mSim.medip.s = post

##################################
#### repeat analysis for medip data after strand assignment using dismiss
i = grep(csScaffold, x = as.character(seqnames(oGRmacs)), perl=T)
oGRdismiss = f_oGRSeparateStrands(oGRmacs[i], csFile = file.choose(), bPaired = T)

## rename sequences
s = gsub("gi\\|\\d+\\|ref\\|(.+)\\|", replacement = '\\1', x = as.character(seqnames(oGRdismiss)), perl=T)
oGRdismiss.s = GRanges(s, ranges(oGRdismiss), strand=strand(oGRdismiss), mcols=mcols(oGRdismiss))

## calculate overlaps
### for each defined genomic feature, count how many peaks overlap
prom = overlapsAny(oGRprom, oGRdismiss.s)
ds = overlapsAny(oGRds, oGRdismiss.s)
genes = overlapsAny(oGRgene, oGRdismiss.s)
exon = overlapsAny(oGRexon, oGRdismiss.s)
intron = overlapsAny(oGRintron, oGRdismiss.s)

mDat = cbind(prom=as.numeric(table(prom)), exon=as.numeric(table(exon)), intron=as.numeric(table(intron)),
             ds=as.numeric(table(ds)))

rownames(mDat) = c('False', 'True')

ivDat = colSums(mDat)
# convert to rates
ivProb.dismiss.s = sweep(mDat, 2, ivDat, '/')['True',] * 1000

#prior = rbeta(1000, ivDat['True'], ivDat['False'])
post = sapply(ivProb.dismiss.s, function(x) {
  r = rgamma(1000, shape = jef.prior['alpha'] + sum(x) , jef.prior['beta']+ length(x))
  return(r)
})

mSim.dismiss.s = post


##### repeat on bisulphite dataset
## choose the correct bismark methyl extractor files
# read the Original Top (OT) file first
f = file.choose()
f2 = file.choose()
oGR.chh.ot = f_oGRReadBismarkMethylExtractor(file = f, '+')
oGR.chh.ob = f_oGRReadBismarkMethylExtractor(file = f2, '-')
f = file.choose()
f2 = file.choose()
oGR.chg.ot = f_oGRReadBismarkMethylExtractor(file = f, '+')
oGR.chg.ob = f_oGRReadBismarkMethylExtractor(file = f2, '-')
f = file.choose()
f2 = file.choose()
oGR.cg.ot = f_oGRReadBismarkMethylExtractor(file = f, '+')
oGR.cg.ob = f_oGRReadBismarkMethylExtractor(file = f2, '-')
# this may give a warning which we can ignore
oGRbis = c(oGR.cg.ob, oGR.cg.ot, oGR.chg.ob, oGR.chg.ot, oGR.chh.ob, oGR.chh.ot)
oGRbis = sort(oGRbis)

## shorten names
s = gsub("gi\\|\\d+\\|ref\\|(.+)\\|", replacement = '\\1', x = as.character(seqnames(oGRbis)), perl=T)
oGRbis.s = GRanges(s, ranges(oGRbis), strand=strand(oGRbis))

# save object and load from here in future
save(oGRbis.s, file='Bee_data/Objects/oGRbis.s')

## if we want to remove low count potentially noisy data
w = countOverlaps(unique(oGRbis.s), oGRbis.s)
f = w <= bs.cut
gr = unique(oGRbis.s)
gr = gr[!f]
f = overlapsAny(oGRbis.s, gr)
oGRbis.s = oGRbis.s[f]
length(oGRbis.s)
w = countOverlaps(unique(oGRbis.s), oGRbis.s)
summary(w)
gc()

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

ivDat = colSums(mDat)
ivProb.bis.s = sweep(mDat, 2, ivDat, '/')['True',] * 1000

#prior = rbeta(1000, ivDat['True'], ivDat['False'])
post = sapply(ivProb.bis.s, function(x) {
  r = rgamma(1000, shape = jef.prior['alpha'] + sum(x) , jef.prior['beta']+ length(x))
  return(r)
})

mSim.bis.s = post

ivProb.bis.s = colMeans(mSim.bis.s)
ivProb.medip.s = colMeans(mSim.medip.s)
ivProb.dismiss.s = colMeans(mSim.dismiss.s)

##### plot bs, medip and dismiss together
c = grey.colors(3)
mBar = rbind(ivProb.medip.s, ivProb.dismiss.s, ivProb.bis.s)
rownames(mBar) = c('MeDIP', 'DISMISS', 'BS-Seq')
l2 = barplot(mBar, beside=T, ylim=c(0, 500), col=c)
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


