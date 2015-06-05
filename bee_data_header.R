# File: bee_data_header.R
# Auth: u.niazi@imperial.ac.uk
# DESC: header file for user defined functions - source before calling other scripts script
# Date: 02/06/2015


##### header files
if (!require(GenomicAlignments) || !require(GenomicFeatures) || !require(rtracklayer)) { 
  stop(paste('Bioconductor libraries GenomicAlignments, GenomicFeatures and rtracklayer required\n
             visit http://www.bioconductor.org/install/ for instructions.'))}

source('dismiss_header.R')

##### PATHS
pWD = getwd()

##### functions used by script

# Function: f_oGRReadBismarkMethylExtractor
# Desc: as input it takes the name of the file (tab separated) created by methyl_extractor script
#       in bismark and the strand i.e. + for OT files (Original Top) and - for OB. returns the data
#       in a GRagnes sorted object
# Args: file: name of bismark file; strand: strand + or -
# Rets: single stranded, sorted GRanges object
f_oGRReadBismarkMethylExtractor = function(file, strand){
  # read the data from the tab separated file
  # also skip the first line as it has comments
  require(GenomicRanges)
  df = read.csv(file, header = F, sep = '\t', skip=1)
  # remove cytosines that are not methylated 
  # they will be marked by -
  df = df[df$V2 != '-',]
  # create ranges of what is left behind
  # it will be plus for original top and minus for original bottom
  gr = GRanges(df$V3, IRanges(df$V4, df$V4), strand=strand)
  gr = sort(gr)
  rm(df)
  # garbage collector
  gc(verbose = FALSE)
  return(gr)  
} # function

# Function: f_oGRRemovePoissonNoiseFromBismark
# Desc: once the Cs that have been methylated and extracted using f_oGRReadBismarkMethylExtractor
#       we have a set of Cs that have been seen as methylated at least once, which may be just noise
#       this function takes that particular range e.g. CHH_OT; assumes that the counts of Cs are poisson
#       and removes any counts under 0.05 p-value under the noise model
# Args: GRanges object created via f_oGRReadBismarkMethylExtractor e.g. CHH_OT
# Rets: GRanges object with possible noise or Cs with low counts removed.
f_oGRRemovePoissonNoiseFromBismark = function(gr){
  # get unique regions 
  gr.unique = unique(gr)
  # count how many events occured over each region
  vec = countOverlaps(gr.unique, gr)
  # take log of data to perform calculations
  vec.log = log(vec) # there should not be any zeroes here
  cut.pt = qpois(0.05, lambda = mean(vec.log), lower.tail = F)
  #cut.pt = round(exp(cut.pt), 0)
  # use this cut.pt to extract ranges from the unique ranges
  f = vec.log < cut.pt
  gr.unique = gr.unique[!f]
  # use these ranges to extract data from original ranges
  f = overlapsAny(gr, gr.unique)
  gr = gr[f]
  return(gr)  
}

# Function: f_LoadObject
# Desc: Takes a file name of R object and returns the object to be assigned to a new variable see
#       http://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
# Args: r.obj.file = r object file name
# Rets: returns the object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}