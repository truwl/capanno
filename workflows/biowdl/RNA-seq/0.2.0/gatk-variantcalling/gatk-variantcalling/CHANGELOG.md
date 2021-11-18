Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 3.0.0
-----------------
+ Updated default images.
+ Estimate the time needed for each scatter by checking the size of the input
  BAM file.
+ Use [scatter-regions](https://github.com/biowdl/chunked-scatter) to replace
  biopet-scatterregions. This allows the pipeline to work with scattersizes
  greater than 2 billion.
+ Added  `bcftools stats` task to generate stats on
  called VCF files.
+ Tasks were updated to contain the `time_minutes` runtime attribute and
  associated `timeMinutes` input, describing the maximum time the task will
  take to run.
+ Refactoring of the pipeline:
    + Split up the pipeline into a single sample variant calling pipeline and 
      a part that performs the joint genotyping. This allows for more elegantly
      integrating the pipeline into other pipelines.
    + Merge steps are only performed when there is more than one scatter. 
      This prevents data from being written twice unnecessarily.
    + `multisample-variantcalling.wdl` is a reference implementation.

version 2.0.0
-----------------
+ Add a scatterSizeMillions parameter to make it easier to set larger scatter 
  sizes.
+ Multisample VCFs are only produced when joint genotyping is used.
+ Add option to output single-sample GVCFs
+ Make Joint Genotyping by GenotypeGVCF an optional step, so the pipeline can 
  be used for RNA variant calling.
+ Make using a dbsnp VCF file optional.
+ Added gender-aware capabilities to the pipeline. This has changed the input
  format.
+ Added inputs overview to the docs.
+ Added parameter_mets.
+ Added wdl-aid to linting.
+ Added miniwdl to linting.

version 1.0.0
---------------------------
+ Combine the bam-to-gvcf and joint-genotyping pipeline into one.
