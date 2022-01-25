---
layout: default
title: Home
---

This workflow can be used to perform expression quantification for multiple
BAM files. Expression levels will be determined for each BAM file/sample
and will be merged together into a single table including all samples.

Expression quantification will be performed using
[StringTie](https://ccb.jhu.edu/software/stringtie/) and
[HTSeq-Count](http://htseq.readthedocs.io/en/master/count.html).

This workflow is part of [BioWDL](https://biowdl.github.io/)
developed by the SASC team
at [Leiden University Medical Center](https://www.lumc.nl/).

## Usage
This workflow can be run using
[Cromwell](http://cromwell.readthedocs.io/en/stable/):

First download the latest version of the workflow wdl file(s)
from the
[github page](https://github.com/biowdl/expression-quantification).

The workflow can then be started with the following command:
```bash
java \
    -jar cromwell-<version>.jar \
    run \
    -o options.json \
    -i inputs.json \
    multi-bam-quantify.wdl
```

Where `options.json` contains the following json:
```json
{
    "final_workflow_outputs_dir": "/path/to/outputs",
    "use_relative_output_paths": true,
    "final_workflow_log_dir": "/path/to/logs/folder"
}
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
    "MultiBamExpressionQuantification.bams": "A list of BAM files and sample identifiers (see 'BAM file input' below).",
    "MultiBamExpressionQuantification.strandedness": "The strandedness of the samples: FR (forward-reverse), RF (reverse-forward) or None.",
    "MultiBamExpressionQuantification.outputDir": "The path to the output directory.",
    "MultiBamExpressionQuantification.referenceGtfFile": "The path to the annotations GTF file. If not specified, Stringtie will be run unguided and the GTF file it produces will be used for HTSeq-Count.",
}
```

#### BAM file input
BAM files need to be given as a list with one item per sample. Each of the
items should be an object containing a `"Left"` element (the sample id) and a
`"Right"` element (the BAM file and its index) following the structure as
shown here:
```json
{
    "left": "Sample identifier",
    "right": {
        "file": "The path to the sample's BAM file",
        "index": "The path to the index for the sample's BAM file"
    }
}
```

#### Example
```json
{
    "MultiBamExpressionQuantification.bams": [
        {
            "left": "s1",
            "right": {
                "file": "/home/user/mapping/results/s1.bam",
                "index": "/home/user/mapping/results/s1.bai"
            }
        },
        {
            "left": "s2",
            "right": {
                "file": "/home/user/mapping/results/s2.bam",
                "index": "/home/user/mapping/results/s2.bai"
            }
        }
    ],
    "MultiBamExpressionQuantification.strandedness": "FR",
    "MultiBamExpressionQuantification.outputDir": "/home/user/expression/results",
    "MultiBamExpressionQuantification.referenceGtfFile": "/home/user/genomes/human/features/ensembl87.gtf"
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
The `multi-bam-quantify` workflow produces two directories:
- **stringtie**: Contains the Stringtie output. Includes two additional files:
  `all_samples.FPKM` and `all_samples.TPM`, which contain the FPKM and
  TPM values for all samples.
- **fragments_per_gene**: Contains the HTSeq-Count output. Also contains a
  file called `all_samples.fragments_per_gene`, which contains the counts
  for all samples.

## Contact
<p>
  <!-- Obscure e-mail address for spammers -->
For any questions about running this workflow and feature requests (such as
adding additional tools and options), please use the
<a href="https://github.com/biowdl/expression-quantification/issues">github issue tracker</a>
or contact the SASC team directly at: 
<a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;">
&#115;&#97;&#115;&#99;&#64;&#108;&#117;&#109;&#99;&#46;&#110;&#108;</a>.
</p>
