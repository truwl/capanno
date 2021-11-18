---
layout: default
title: Home
---

This repository contains the [BioWDL](https://github.com/biowdl)
workflow which can be used for quality control preprocessing and 
reporting of sequencing data.

These workflows are part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team at [Leiden University Medical Center](https://www.lumc.nl/).

## Usage

`QC.wdl` can be run using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):
```
java -jar cromwell-<version>.jar run -i inputs.json QC.wdl
```

### Input

Inputs are provided through a JSON file. The minimally required inputs are
described below, but additional inputs are available.
A template containing all possible inputs can be generated using
Womtool as described in the
[WOMtool documentation](http://cromwell.readthedocs.io/en/stable/WOMtool/).
For an overview of all available inputs, see [this page](./inputs.html).

```JSON
{
  "QC.read1": "Path to file with forward reads / unpaired reads"
}
```
`QC.read1`  is the only required input. In case of read pairs the reverse
read can be set with `QC.read2`. 

Optional inputs:
```JSON
{
  "QC.read2": "Path to file with reverse reads",
  "QC.adapterForward":  "The adapter for the forward reads (read1), default = \"AGATCGGAAGAG\"",
  "QC.adapterReverse": "The adapter for the reverse reads (read2), default = \"AGATCGGAAGAG\")",
  "QC.contaminations": "A list of contaminations to be cleaned with cutadapt (Optional)",
  "QC.runAdapterClipping": "Can be set to false to prevent cutadapt from running.",
  "QC.readgroupName": "What basename should be used to save the fastq files. By default will use the name of the fastq as in <name>.fq.gz",
}
```

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
Alternatively an output directory can be set with `QC.outputDir`.
`QC.outputDir` must be mounted in the docker container. Cromwell will
need a custom configuration to allow this.

#### Example

An example of an inputs.json might look like this:
```JSON
{
  "QC.read1":"/home/user/samples/sample_1/lib_1/rg_1/R1.fq.gz",
  "QC.read2":"/home/user/samples/sample_1/lib_1/rg_1/R2.fq.gz",
  "QC.adapterForward": ["AGATCGGAAGAG"],
  "QC.adapterReverse": ["AGATCGGAAGAG"]
}
```

Note that `adapterBoth` uses a list of strings instead of a single string.
This is because cutadapt accepts multiple adapters.

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

A new set of FASTQ files from which detected adapters have been clipped and a
set of quality reports.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any question related to these workflows, please use the
<a href='https://github.com/biowdl/QC/issues'>github issue tracker</a>
or contact the SASC team directly at: 
<a href='&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;'>
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
