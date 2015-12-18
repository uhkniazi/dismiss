# dismiss
DISMISS is an R script, which as an additional step in MeDIP-Seq data analysis workflow, enables the allocation of strands to methylated DNA regions. It does this by analyzing the proportions of first mate reads aligning to the methylated locus from the plus and minus strands.

##Requirements

Dismiss is a collection of functions implemented in R  and require the following to be present on the system:

<ul>
<li>I.	R (version 3.1.2)</li>
<li>II.	Bioconductor ,  (version 3.0) with packages</li>
  <ul>
    <li>A.  GenomicAlignments</li>  
    <li>B.	GenomicFeatures</li>  
    <li>C.	rtracklayer</li>
  </ul>
</ul>

##Running dismiss

To run dismiss, the following files should be in the current path or directory:
    <ul>
      <li>`dismiss_header.R` (the script with all the required functions)</li>
      <li>`dismiss_macs2_extractor.R` (utility script to process MACS2 results)</li>
      <li>Tab separated `xls` results file for methylated regions for peak calling performed by MACS2. </li>
      <li>Alignment file in `.bam` format (the same file used to perform peak calling)</li>
      <li>Bam index file with `.bai` extension (if bam file is called `file.bam` then index file should be called `file.bam.bai`)</li>
  </ul>
  
From the GNU/Linux command line the following command will run dismiss assuming that the files above are in your current working directory;

```
Rscript dismiss_macs2_extractor.R MACS2_result.xls File.bam [p or s] 
```
  
where the last argument is essential and should be set to either;  
    p = for paired end data OR  
    s = for single end data  

##Additional Information
  
<ul>
    <li>dismiss on Galaxy http://share-galaxy.ibers.aber.ac.uk/dismiss</li>
    <li>Manual http://uhkniazi.github.io/dismiss</li>
</ul>
