---
layout: default
title: Home
---

This workflow can be used to collect a variety of metrics from a BAM file.
Metrics are collected using picard and samtools.

This workflow is part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team
at [Leiden University Medical Center](https://www.lumc.nl/).

## Usage
This workflow can be run using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):

First download the latest version of the workflow wdl file(s)
from the
[github page](https://github.com/biowdl/BamMetrics).


The workflow can then be started with the following command:
```bash
java \
    -jar cromwell-<version>.jar \
    run \
    -o options.json \
    -i inputs.json \
    bammetrics.wdl
```

### Inputs
Inputs are provided through a JSON file. The minimally required inputs are
described below, but additional inputs are available.
A template containing all possible inputs can be generated using
Womtool as described in the
[WOMtool documentation](http://cromwell.readthedocs.io/en/stable/WOMtool/).
For an overview of all available inputs, see [this page](./inputs.html).

```json
{
    "BamMetrics.referenceFasta": "The path to the reference fasta file.",
    "BamMetrics.referenceFastaFai": "The path to the index for the reference fasta.",
    "BamMetrics.referenceFastaDict": "The path to the sequence dictionary dict file for the reference fasta.",
    "BamMetrics.bam": "A path to an input BAM file.",
    "BamMetrics.bamIndex": "A path to the index of the BAM file."
}
```

In the case of RNAseq samples the following inputs are also required. By
providing these inputs additional RNAseq metrics will be collected.

```json
{
    "BamMetrics.refRefflat": "A path to a Refflat annotation file for the reference.",
    "BamMetrics.strandedness": "The strandedness of an RNAseq sample. This should be on of 'None', 'FR' (forward-reverse) or 'RF' (reverse-forward), defaults to 'None'."
}
```

For targeted sequencing samples the following inputs are also required. By
providing these inputs additional targeted PCR metrics will be collected.

```json
{
    "BamMetrics.targetIntervals": "A list of paths to the bed files containing the target regions.",
    "BamMetrics.ampliconIntervals": "The path to the bed file containing the amplicon regions."
}
```

An output directory can be set using an `options.json` file. See the
[cromwell documentation](
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

Alternatively an output directory can be set with `BamMetrics.outputDir`.
`BamMetrics.outputDir` must be mounted in the docker container. Cromwell will
need a custom configuration to allow this.

#### Example
The following is an example of what an inputs JSON might look like:
```json
{
    "BamMetrics.referenceFasta": "/home/user/genomes/human/GRCh38.fasta",
    "BamMetrics.referenceFastaFai": "/home/user/genomes/human/GRCh38.fasta.fai",
    "BamMetrics.referenceFastaDict": "/home/user/genomes/human/GRCh38.dict",
    "BamMetrics.outputDir": "/home/user/analysis/metrics",
    "BamMetrics.bam": "/home/user/mapping/results/s1.bam",
    "BamMetrics.bamIndex": "/home/user/mapping/results/s1.bai",
    "BamMetrics.targetIntervals": [
        "/home/user/analysis/target1.bed",
        "/home/user/analysis/target2.bed"
    ],
    "BamMetrics.ampliconIntervals": "/home/user/analysis/amplicon.bed"
}
```

## Dependency requirements and tool versions
Biowdl workflows use docker images to ensure reproducibility. This
means that biowdl workflows will run on any system that has docker
installed. Alternatively they can be run with singularity.

For more advanced configuration of docker or singularity please check
the [cromwell documentation on containers](
https://cromwell.readthedocs.io/en/stable/tutorials/Containers/).

Images from [biocontainers](https://biocontainers.pro) are preferred for
biowdl workflows. The list of default images for this workflow can be
found in the default for the `dockerImages` input.

## Output
A directory containing the various metrics collected.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any questions about running this workflow and feature requests (such as
adding additional tools and options), please use the
<a href="https://github.com/biowdl/bammetrics/issues">github issue tracker</a>
or contact the SASC team directly at: 
<a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;">
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
