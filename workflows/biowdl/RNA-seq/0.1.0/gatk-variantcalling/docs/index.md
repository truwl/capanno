---
layout: default
title: Home
---

This workflow can be used to generate a multisample VCF file from BAM 
files using GATK HaplotypeCaller.

This workflow is part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team at [Leiden University Medical Center](
https://www.lumc.nl/).

## Usage
This workflow can be run using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):
```
java -jar cromwell-<version>.jar run -i inputs.json multisample-variantcalling.wdl
```

The pipeline can be integrated into other pipelines as well. It has been split
up into three parts for convenience.

+ `calculate-regions.wdl`: Calculates the scatters taking into account 
  (optional) the non-PAR regions and regions of interest.
+ `single-sample-variantcalling.wdl`: Variant calling using the given scatters
  and non-PAR regions. If no scattering is used and there are no non-PAR 
  regions given, this is the same as a plain call to HaplotypeCaller. Can be
  run in GVCF mode.
+ `jointgenotyping.wdl`. Runs the jointgenotyping (scattered). This step is
  optional. Requires the single sample pipeline to be run in GVCF mode.

### Inputs
Inputs are provided through a JSON file. The minimally required inputs are
described below and a template containing all possible inputs can be generated
using Womtool as described in the
[WOMtool documentation](http://cromwell.readthedocs.io/en/stable/WOMtool/).
For an overview of all available inputs, see [this page](./inputs.html).
```json
{
  "MultisampleCalling.referenceFasta": "A reference fasta file",
  "MultisampleCalling.referenceFastaFai": "The index for the reference fasta",
  "MultisampleCalling.referenceFastaDict": "The dict file for the reference fasta",
  "MultisampleCalling.bamFilesAndGenders": "A list of structs. Each struct contains the bam file, the index and the gender of the sample. The gender is optional. " 
 }
```
`gender` can be set to `null`. Actionable values are `female`, `male`, `F`, 
`f`, `M`, or `m`. Any other values and `null` will handle the sample as an 
"unknown" gender.

The following actions are taken for each gender:

+ `female`: The non-PAR regions of X will be called with a ploidy of 2. Y will 
  not be called.
+ `male`: The non-PAR regions of both X and Y will be called with a ploidy of 1. 
+ `unknown`: The non-PAR regions of X will be called with a ploidy of 2. The 
   non-PAR region of Y will be called with a ploidy of 1.

Some additional inputs which may be of interest are:
```json
{
  "MultisampleCalling.dbsnpVCF": "A dbSNP VCF file with known variants.",
  "MultisampleCalling.dbsnpVCFIndex": "Index (.tbi) for the dbSNP VCF file",
  "MultisampleCalling.regions": "The path to a bed file containing the regions for which variant calling will be performed",
  "MultisampleCalling.scatterSize": "The size of scatter regions (see explanation of scattering below), defaults to 10,000,000",
  "MultisampleCalling.vcfBasename": "The basename of the to be outputed VCF files, defaults to 'multisample'",
  "MultisampleCalling.XNonParRergions": "Bed file with the non-PAR regions of X. Required for gender-aware variant calling.",
  "MultisampleCalling.YNonParRegions": "Bed file with the non-PAR regions of Y. Required for gender-aware variant calling.",
  "MultisampleCalling.singleSampleGvcf": "Output Gvcfs for every single sample."
}
```
When the X and Y non-PAR regions are not both provided, GATK will call all 
chromosomes with ploidy 2 naively. 

By default GVCFs are not created for every single sample as it requires a lot
of extra space and write actions. If the option to output single sample GVCFs
is turned on, this will happen in parallel with the creation of the multisample
GVCF file.

An output directory can be set using an `options.json` file. See [the
cromwell documentation](
https://cromwell.readthedocs.io/en/stable/wf_options/Overview/) for more
information.

Example `options.json` file:
```JSON
{
"final_workflow_outputs_dir": "my-analysis-output",
"use_relative_output_paths": true,
"default_runtime_attributes": {
  "docker_user": "$EUID"
  }
}
```
Alternatively an output directory can be set with `MultisampleCalling.outputDir`.
`MultisampleCalling.outputDir` must be mounted in the docker container. Cromwell will
need a custom configuration to allow this.

#### Example
```json
{
  "MultisampleCalling.dbsnpVCF": "/home/user/genomes/human/dbsnp/dbsnp-151.vcf.gz",
  "MultisampleCalling.dbsnpVCFIndex": "/home/user/genomes/human/dbsnp/dbsnp-151.vcf.gz.tbi",
  "MultisampleCalling.referenceFasta": "/home/user/genomes/human/GRCh38.fasta",
  "MultisampleCalling.referenceFastaFai": "/home/user/genomes/human/GRCh38.fasta.fai",
  "MultisampleCalling.referenceFastaDict": "/home/user/genomes/human/GRCh38.dict",
  "MultisampleCalling.vcfBasename": "s1",
  "MultisampleCalling.XNonParRegions": "/home/user/genomes/human/x_non_par.bed",
  "MultisampleCalling.YNonParRegions": "/home/user/genomes/human/y_non_par.bed",
  "MultisampleCalling.outputDir": "/home/user/analysis/results/",
  "MultisampleCalling.bamFilesAndGenders": [
    {"file": "/home/user/mapping/results/s1_1.bam",
     "index":  "/home/user/mapping/results/s1_1.bai",
     "gender":"male"},
    {"file": "/home/user/mapping/results/s1_2.bam",
     "index": "/home/user/mapping/results/s1_2.bai",
     "gender": "female"},
    {"file": "/home/user/mapping/results/s1_3.bam",
     "index":  "/home/user/mapping/results/s1_3.bai",
     "gender": null}
   ]
}
```

### Dependency requirements and tool versions
Biowdl pipelines use docker images to ensure  reproducibility. This
means that biowdl pipelines will run on any system that has docker
installed. Alternatively they can be run with singularity.

For more advanced configuration of docker or singularity please check
the [cromwell documentation on containers](
https://cromwell.readthedocs.io/en/stable/tutorials/Containers/).

Images from [biocontainers](https://biocontainers.pro) are preferred for
biowdl pipelines. The list of default images for this pipeline can be
found in the default for the `dockerImages` input.

### output
A multisample vcf file and a multisample gvcf file.

## Scattering
This pipeline performs scattering to speed up analysis on grid computing
clusters. This is done by splitting the reference genome into regions of
roughly equal size (see the `scatterSize` input). Each of these regions will be
analyzed in separate jobs, allowing them to be processed in parallel.

For each BAM file input in the pipeline these scatters will be used. For 
example, with 5 BAM files and 4 scatters, 20 jobs will run to genotype the BAM
files.
The resulting GVCF files will be merged. Then GatkGenotypeGVCF will be 
scattered over the 4 scatters to create a VCF per scatters.
The output is merged into a multisample vcf.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any question related to this workflow, please use the
<a href='https://github.com/biowdl/gatk-variant-calling/issues'>github issue tracker</a>
or contact the SASC team directly at: 
<a href='&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;'>
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
