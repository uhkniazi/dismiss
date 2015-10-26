# Copyright (C) 2014-2015  Umar Niazi
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# File: dismiss_header.R
# Auth: umn@aber.ac.uk, uhkniazi@gmail.com
# DESC: header file for user defined functions - source before calling other scripts script
# Date: 24/10/2014
# NOTE: call function f_oGRSeparateStrands if using your own main scripts for strand separation

# load the libraries 
if (!require(GenomicAlignments) || !require(GenomicFeatures) || !require(rtracklayer)) { 
  stop(paste('Bioconductor libraries GenomicAlignments, GenomicFeatures and rtracklayer required\n
             visit http://www.bioconductor.org/install/ for instructions.'))}


## functions used by script

# NAME: f_bComparePoissonLiklihoodApply
# ARGS: a vector with 2 components
# DESC: if the event i.e. total number of reads on a signal interval can be divided into
#       2 subevents - then we can see if the difference between these 2 subevents 
#       is significant or not
# RETS: TRUE - if the 2 rates are equal
f_bComparePoissonLiklihoodApply = function(ivVector, cutoff=0.2) {
  # check page 293 of introduction to statistical thought
  # for details of how to divide a poisson event into 2
  iLam.1 = ivVector[1]
  iLam.2 = ivVector[2]
  iLam.total = iLam.1 + iLam.2
  # the event lambda total is composed of 2 subevents
  # these subevents have a binomial distribution
  x = 0:iLam.total
  # liklihood of first distribution
  lik.1 = dpois(iLam.1, lambda=x)
  lik.1 = lik.1/max(lik.1)
  lik.2 = dpois(iLam.2, lambda=x)
  lik.2 = lik.2/max(lik.2)
  # if any of the values match or intersect between the 2 
  # distributions then we can say they are not different from
  # one another
  x.1 = x[which(lik.1 > cutoff)]
  x.2 = x[which(lik.2 > cutoff)]
  if (length(x.1[x.1 %in% x.2]) > 0) return(TRUE) else return(FALSE)
} #function

# Function: f_mFastCountMateReadsOverRanges
# Desc: it is a fast version of the function CountMateReadsOverRanges, but will use more memory
#       it uses the granges object (typically calculated by a peak caller like MACS)
#       and a bam file name, and counts how many reads aligned from first mate and both mates
#       over those ranges.
# Args: gr.signal: granges object containing the coordinates for signals in genome
#       bam.file: the bam file used for alignment and peak calling
#       bPaired: boolean value, default is true. if data was single end read then set it false
#                as in single end reads there are no first mate reads
# Rets: a matrix N X 4 size (N being length of GRanges signal object)
f_mFastCountMateReadsOverRanges = function(gr.signal, bam.file, bPaired=T){
  # matrix to store return object
  mRets = matrix(NA, nrow=length(gr.signal), ncol=4, dimnames=list(NULL, c('mp', 'mm', 'bp', 'bm')))
  # set flag to read bam file first mates
  flag = scanBamFlag(isFirstMateRead=T)
  if (bPaired){  bam = readGAlignments(bam.file, param=ScanBamParam(flag=flag, which=gr.signal))
                 # how many reads align with + strand and - strands
                 mRets[,'mp'] = countOverlaps(gr.signal, bam[strand(bam)=='+'])
                 mRets[,'mm'] = countOverlaps(gr.signal, bam[strand(bam)=='-'])
  } #if
  # now read data coming from both mates 
  bam = readGAlignments(bam.file, param=ScanBamParam(which=gr.signal))
  # how many reads align with + strand and - strands
  mRets[,'bp'] = countOverlaps(gr.signal, bam[strand(bam)=='+'])
  mRets[,'bm'] = countOverlaps(gr.signal, bam[strand(bam)=='-'])
  # if data was single end read then first mates will be zero in matrix
  if (bPaired == F){ mRets[,'mp'] = mRets[,'bp']
                     mRets[,'mm'] = mRets[,'bm'] } # if
  return(mRets)  
} # function


# Function: f_dfReadMACS2xls
# Desc: Simple function to read data from MACS2 xls file output that is tab separated
# Args: xls.file: path to xls file produced by MACS2
#       skip.lines: by default MACS2 has 22 lines of comments
#       separator: by default it is a tab separated file
# Rets: a data frame with data
f_dfReadMACS2xls = function(xls.file, separator='\t'){
  # read the data in MACS2 xls file
  # it is a tab separated file with n lines of comments
  # the comment lines start with a #, so skip these lines
  infile = file(xls.file, 'rt')
  input = readLines(infile, n = 1000)
  n = grep('^#', x = input)
  close(infile)
  # we have n lines of comments, so skip n+1 lines  
  df = read.csv(file=xls.file, header=T, sep=separator, skip=length(n)+1)
  return(df)
} # function


# Function: f_oGRMACS2toGRanges
# Desc: Simple function to read data from MACS2 xls file output that is tab separated
# Args: xls.file: path to xls file produced by MACS2
#       skip.lines: by default MACS2 has 22 lines of comments
#       separator: by default it is a tab separated file
# Rets: a Granges object of MACS2 xls data
f_oGRMACS2toGRanges = function(xls.file, separator='\t'){
  # use previously defined function to read data frame
  df = f_dfReadMACS2xls(xls.file, separator)
  # convert data frame to GRanges object
  gr = GRanges(seqnames=as.character(df$chr), ranges=IRanges(df$start, df$end))
  dt = DataFrame(abs_summit=df$abs_summit, pileup=df$pileup, 
                 fold_enrichment=df$fold_enrichment, neg.log10.qval=df$X.log10.qvalue.)
  mcols(gr) = dt
  return(gr)
} # function

# Function: f_WriteGff3
# Desc: Takes a GRanges object, separates it into 3 groups based on strands (*, + or -)
#       checks for length of each group, if length is greater than 0
#       writes that object into a gff3 format using rtracklayer library
# Args: a granges object to write
# Rets: none
f_WriteGff3 = function(oGR){
  library(rtracklayer)
  # write 3 gff files for *, + and - strand
  if (length(oGR[strand(oGR)=='*']) > 0){
    name.out = paste('oGRdismiss_', make.names(date()), '_double_st.gff3', sep='')
    export(oGR[strand(oGR)=='*'], name.out, format='gff3')
  }  
  if (length(oGR[strand(oGR)=='+']) > 0){
    name.out = paste('oGRdismiss_', make.names(date()), '_plus_st.gff3', sep='')
    export(oGR[strand(oGR)=='+'], name.out)
  }  
  if (length(oGR[strand(oGR)=='-']) > 0){
    name.out = paste('oGRdismiss_', make.names(date()), '_minus_st.gff3', sep='')
    export(oGR[strand(oGR)=='-'], name.out, format='gff3')
  }
} # function end 


# Function: f_oGRSeparateStrands
# Desc: main function to separate strands, makes use of other previous functions
# Args: oGRsignal.strand = a granges object with ranges for locations of peaks
#       csFile = location of bam file, it assumes a bam index file is present
#       bPaired = boolean value, set to true for paired end reads or else set to false 
# Rets: a Granges object of MACS2 xls data
f_oGRSeparateStrands = function(oGRsignal.strand, csFile, bPaired){
  # count reads aligning to regions
  mFirstMate = f_mFastCountMateReadsOverRanges(oGRsignal.strand, bam.file = csFile, bPaired)
  mp = mFirstMate[,'mp'] + 1 # add a one to avoid zero probabilities and divisions by zero
  mm = mFirstMate[,'mm'] + 1
  n = mp + mm
  p.plus = mp/n
  p.minus = 1-p.plus
  mFirstMate = cbind(mFirstMate, p.plus, p.minus)
  ## assign strands to regions
  # check if the 2 reads on plus and minus strands are equal
  # returns a boolean true if equal else false
  t = apply(mFirstMate[,c('mp', 'mm')], 1, f_bComparePoissonLiklihoodApply)
  # those regions where t is true are double stranded
  # others will be plus stranded if reads on plus side > reads on minus side
  t.f = rep(NA, length.out=nrow(mFirstMate))
  for (i in 1:length(t)){
    if (t[i]){
      t.f[i] = 2 # both strands
    } else if (mFirstMate[i,'mp'] > mFirstMate[i,'mm']) {
      t.f[i] = 1 # else if plus strand gives more first mates
    } else t.f[i] = 0 # else minus strand gives more first mates
  } #for
  st = factor(t.f, labels=c('-', '+', '*'))
  # add data to GRanges object
  ## add data to GRanges object
  oGRsignal.strand$strand_fac = t.f
  strand(oGRsignal.strand) = st
  df = mcols(oGRsignal.strand)
  df = cbind(df, DataFrame(mFirstMate))
  mcols(oGRsignal.strand) = df
  return(oGRsignal.strand)  
}