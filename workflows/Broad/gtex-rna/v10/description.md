
This repository contains all components of the RNA-seq pipeline used by the GTEx Consortium, including alignment,
expression quantification, and quality control.

- SamToFastq: BAM to FASTQ conversion
- STAR: spliced alignment of RNA sequence reads (v2.5.3a)
- Picard MarkDuplicates: mark duplicate reads
- RSEM transcript expression quantification (v1.3.0)
- bamsync: utility for transferring QC flags from the input BAM and for re-generating read group IDs
- RNA-SeQC: QC metrics and gene-level expression quantification (v1.1.9)

Reference indexes for STAR and RSEM are needed to run the pipeline. All reference files are available at [gs://gtex-resources](https://console.cloud.google.com/storage/browser/gtex-resources).

GTEx releases from V8 onward are based on the GRCh38/hg38 reference genome. Please see [TOPMed_RNAseq_pipeline.md](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) for details and links for this reference. Releases up to V7 were based on the GRCh37/hg19 reference genome ([download](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta)).

Release V8 uses the [GENCODE v26](https://www.gencodegenes.org/human/release_26.html) annotation. Releases V6/V6p and V7 used [GENCODE v19](https://www.gencodegenes.org/human/release_19.html).

For hg19-based analyses, the GENCODE annotation should be patched to use Ensembl chromosome names:
```
zcat gencode.v19.annotation.gtf.gz | \
  sed 's/chrM/chrMT/;s/chr//' > gencode.v19.annotation.patched_contigs.gtf
```
A 2x76 bp paired-end sequencing protocol, will use a sjdbOverhang of 75