---
layout: default
title: Home
---

This workflow performs preprocessing steps required for variantcalling based
on the
[GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/).
This workflow can be used for both DNA data and RNA-seq data. It recalibrates
a BAM file and optionally splits spliced reads.

This workflow is part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team at [Leiden University Medical Center](https://www.lumc.nl/).

## Usage
This workflow can be run using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell-<version>.jar run -i inputs.json gatk-preprocess.wdl
```

### Inputs
Inputs are provided through a JSON file. The minimally required inputs are
described below and a template containing all possible inputs can be generated
using Womtool as described in the
[WOMtool documentation](http://cromwell.readthedocs.io/en/stable/WOMtool/).
For an overview of all available inputs, see [this page](./inputs.html).
```json
{
  "GatkPreprocess.referenceFasta": "The path to the reference fasta file",
  "GatkPreprocess.referenceFastaFai": "The path to the index for the reference fasta",
  "GatkPreprocess.referenceFastaDict": "The path to the sequence dictionary dict file for the reference fasta",
  "GatkPreprocess.bam": "A path to an input BAM file",
  "GatkPreprocess.bamIndex": "A path to the index of the BAM file.",
  "GatkPreprocess.bamName": "The name for the output bam. The final output will be <bamName>.bam or <bamName>.bqsr",
  "GatkPreprocess.dbsnpVCF": "A path to a dbSNP VCF file",
  "GatkPreprocess.dbsnpVCFIndex": "The path to the index (.tbi) file associated with the dbSNP VCF"
}
```

Some additional inputs that may be of interest are:
```json
{
  "GatkPreprocess.scatters": "A list of bed files describing the regions to be processed.",
  "GatkPreprocess.splitSplicedReads": "Whether or not SplitNCigarReads should be executed (recommended for RNA-seq data), defaults to false",
}
```
Each bed file supplied with `scatters` will be used in a seperate job for
most of the steps taken in this workflow. This will allow for parallelization
if the backend used supports this. It is recommended to use this input and supply
one bed file per chromosome (small chromosomes can be together in one bed file).

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
Alternatively an output directory can be set with `GatkPreprocess.outputDir`.
`GatkPreprocess.outputDir` must be mounted in the docker container. Cromwell will
need a custom configuration to allow this.

#### Example
```json
{
  "GatkPreprocess.referenceFasta": "/home/user/genomes/human/GRCh38.fasta",
  "GatkPreprocess.referenceFastaFai": "/home/user/genomes/human/GRCh38.fasta.fai",
  "GatkPreprocess.referenceFastadict": "/home/user/genomes/human/GRCh38.dict",
  "GatkPreprocess.bamName": "s1_preprocessed",
  "GatkPreprocess.dbsnpVCF": "/home/user/genomes/human/dbsnp/dbsnp-151.vcf.gz",
  "GatkPreprocess.dbsnpVCFIndex": "/home/user/genomes/human/dbsnp/dbsnp-151.vcf.gz.tbi",
  "GatkPreprocess.bam": "home/user/mapping/results/s1.bam",
  "GatkPreprocess.bamIndex":"/home/user/mapping/results/s1.bai",
  "GatkPreprocess.splitSplicedReads": true
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

### Output
This workflow will produce a BQSR report named according to the `bamName`
input (bamName + '.bqsr'). If one of the `splitSplicedReads` or
`outputRecalibratedBam` inputs is set to true, a new BAM file (bamName +
'.bam') will be produced as well.

## Scattering
This pipeline performs scattering to speed up analysis on grid computing
clusters. This is done by splitting the reference genome into regions of
roughly equal size (see the `scatterSize` input). Each of these regions will
be analyzed in separate jobs, allowing them to be processed in parallel.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any question about running this workflow and feature requests, please use
the
<a href='https://github.com/biowdl/gatk-preprocess/issues'>github issue tracker</a>
or contact
the SASC team
 directly at: 
<a href='&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;'>
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
