This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq and DNase-seq data.
The pipeline can be run on compute clusters with job submission engines as well as on stand alone machines.
It inherently makes uses of parallelized/distributed computing.
Pipeline installation is also easy as most dependencies are automatically installed.
The pipeline can be run end-to-end, starting from raw FASTQ files all the way to peak calling and signal track generation using a single caper submit command.
One can also start the pipeline from intermediate stages (for example, using alignment files as input).
The pipeline supports both single-end and paired-end data as well as replicated or non-replicated datasets. The outputs produced by the pipeline include

1. formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data,
2. analysis of reproducibility
3. stringent and relaxed thresholding of peaks,
4. fold-enrichment and pvalue signal tracks.

The pipeline also supports detailed error reporting and allows for easy resumption of interrupted runs. 
It has been tested on some human, mouse and yeast ATAC-seq datasets as well as on human and mouse DNase-seq datasets.

The ATACseq pipeline protocol specification is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). 
Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

---
####Job Costs*


Examples:

|  Public example     |   Job cost  | Notes |
| ------------------ | ---------------- |-------------------| 
| [PVTT_ATAC_seq](https://truwl.com/workflows/library/ENCODE%20ATAC-seq%20pipeline/v1.8.0/instances/WF_e85df4.f10.9038) | $19.43  | Paired end 2-reps.|
| [ENCSR113MBR ATAC-seq on human adrenal gland](https://truwl.com/workflows/library/ENCODE%20ATAC-seq%20pipeline/v1.8.0/instances/WF_e85df4.f10.e55f) | $23.36  | Paired-end 2-reps. |
| [ENCODE ATAC-seq analysis ENCSR159GFS](https://truwl.com/workflows/library/ENCODE%20ATAC-seq%20pipeline/v1.8.0/instances/WF_e85df4.f10.5238) | $8.72 | Paired-end 1-rep.|

*Job cost examples are for estimates only. To get a more accurate idea of job costs, try running a single job before running many jobs.