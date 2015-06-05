# dismiss
DISMISS is an R script, which as an additional step in MeDIP-Seq data analysis workflow, enables the allocation of strands to methylated DNA regions. It does this by analyzing the proportions of first mate reads aligning to the methylated locus from the plus and minus strands.

# Bee_data branch
Generate figure for comparison of Bee MeDIP-Seq, DISMISS and BS-Seq data

# bee_data_header.R
header file for user defined functions - source before calling other scripts script. 

# bee_dismiss_vs_bs.R
The script does the comparison between the 3 data sets. It loads the GFF file for the bee ref_Amel_4.5_top_level.gff3 (loaded from
NCBI). Selects the genes and exons and creates granges objects. Reads the MACS2 data using the custom function in dismiss header. As 
the scaffold names in the gff are shorter and they are longer in the fasta files used for alignments, we use gsub to shorten the
names. GRanges features are created i.e. 2k upstream, exons, introns, and 2k downstream. The Dismiss data is obtained by using the
dismiss functions on MACS2 data. The BS-Seq data is imported by a custom function to import bismark data and the low count 5mCs are
removed (i.e. potentially noisy data). The overlaps of features (query) vs methylation data (subject) is done and the conditional 
probability i.e. P(y | Theta) is calculated as a multinomial distribution. The confidence intervals are calculated by simulation. 
We summarize the data as bar plots with error bars for 95% confidence intervals.
