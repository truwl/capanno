This app processes whole genome sequencing (WGS) or whole exome sequencing (WES) data - it aligns raw sequencing data to the human reference genome and additionally performs duplicate marking. 
It generates an OQFE CRAM file from a FASTQ or mapped CRAM file, following the Original Quality Functionally Equivalent (OQFE) protocol, which is a revision of the Functionally Equivalent (FE) protocol. 
The app follows the FE methodology by the adoption of a standard GRCh38 reference genome with alternate loci, read alignment with BWA-MEM, inclusion of supplementary alignments, duplicate marking that includes the supplementary alignments, CRAM compression and restricted tag usage.
  
In accordance with the OQFE protocol, the app retains the original quality scores of the reads, which allows for recovery of original FASTQs from the generated OQFE CRAM, and implements updated versions of the used programs.

---
####Job Costs*

Not enough data to determine cost estimates for this workflow yet. Shouldn't be more than a dollar or two.

*Job cost examples are for estimates only. To get a more accurate idea of job costs, try running a single job before running many jobs.