# dismiss
DISMISS is an R script, which as an additional step in MeDIP-Seq data analysis workflow, enables the allocation of strands to methylated DNA regions. It does this by analyzing the proportions of first mate reads aligning to the methylated locus from the plus and minus strands.

# Bee_data branch
Generate figure for comparison of Bee MeDIP-Seq, DISMISS and BS-Seq data

# bee_data_header.R
header file for user defined functions - source before calling other scripts script. 

# bee_dismiss_vs_bs.R (Redundant) 
## replaced by bee_dismiss_medip_bs_rates.R
The script does the comparison between the 3 data sets. It loads the GFF file for the bee ref_Amel_4.5_top_level.gff3 (loaded from
NCBI). Selects the genes and exons and creates granges objects. Reads the MACS2 data using the custom function in dismiss header. As 
the scaffold names in the gff are shorter and they are longer in the fasta files used for alignments, we use gsub to shorten the
names. GRanges features are created i.e. 2k upstream, exons, introns, and 2k downstream. The Dismiss data is obtained by using the
dismiss functions on MACS2 data. The BS-Seq data is imported by a custom function to import bismark data and the low count 5mCs are
removed (i.e. potentially noisy data). The overlaps of features (query) vs methylation data (subject) is done and the conditional 
probability i.e. P(y | Theta) is calculated as a multinomial distribution. The confidence intervals are calculated by simulation. 
We summarize the data as bar plots with error bars for 95% confidence intervals.

# bee_dismiss_medip_bs_rates.R (Figure 4)
Most of the script is similar to the previous redundant script. The rate of methylation signal overlap per 1000 genomic features of interest is calculated. This was done beacuse the multinomial model did not consider the abundance of features like introns and exons and hence was perhaps not the best way to show the data. The new model can either use a binomial model TRUE/FALSE for each feature category, however some features have a very small proportion of TRUE, and hence the binomial parameter theta can be very small. In such cases a poisson rate variable is more appropriate - and we use a jefferey's non informative prior for the gamma distribution, to calculate the posterior gamma (rate) parameter for each feature, that represents the rate of methylation signal per 1000 features. 