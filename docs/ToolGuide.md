# Introduction

This document contains Truwl's best practices for the describing command
line tools with metadata files and [Common Workflow
Language](https://github.com/common-workflow-language/common-workflow-language)
(CWL) files. Several practices here were adapted from [CommonsWorkflow
Language User Guide: Recommended
Practices](http://www.commonwl.org/user_guide/rec-practices/). Best
practices for describing workflow language files from other workflow
languages will be added as developed.

## Definitions

The words "MUST", "REQUIRED", "SHOULD" and "MAY" have the meanings
described in [Key words for use in RFCs to Indicate Requirement
Levels](https://www.ietf.org/rfc/rfc2119.txt)

## Dividing tools into subtools

All tools in `capanno` are divided into subtools. Subtools are the
subcommands or modes that are specified when the tool is called
e.g. `<tool_name> <subtool_name> <arguments>`. For tools that can be
called without a subcommand, the subtool name is `__main__` and the
files associated with the main tool are placed in a directory with the
same named after the tool name, e.g. of for cat `tools/cat/8.x/cat`.

### Metadata files

There are 2 types of metadata files for tools:

-   **Common metadata**: Metadata that is common to all subtools. The
    main source of metadata for all subtools is typically from common
    metadata.

-   **Subtool metadata**: Metadata in subtool metadata files can add
    metadata specific to the subtool or override metadata from the
    common metadata.

Metadata files can be initialized (along with the proper directory
structure) programmatically as described in [Getting
Started](../Getting_Started.md).

### Scope of tool workflow language files

Workflow language files must be divided into the same subtool structure
as metadata, i.e., there must be only one workflow language file per
subtool, per workflow language, i.e. a single subtool may have one of
each CWL, WDL, Snakemake, and Nextflow files. Subtools typically contain
mutually exclusive arguments that should be divided into separate
workflow language files. For example md5sum CWL file would be divided
into `md5sum/<versionName>/md5sum/md5sum.cwl` and
`md5sum/<versionName>/md5sum_check/md5sum-check.cwl` since
`md5sum --check [File]` has arguments that are not relevant when running
`md5sum` without the --check subcommand.

# Tool metadata files

This describes how to document metadata for command line tools described
in this repository. Metadata files are specified in YAML and can be
generated programmatically (see [getting
started](../../docs/Getting_Started.md)). Most metadata keys are defined
by [schema.org](https://schema.org/) vocabularies and can be adapted to
be imported by CWL files.

## List order

Metadata fields that accept lists as values should be ordered by
importance if applicable.

## Metadata sources

We recommend populating metadata files from
[bio.tools](https://bio.tools/) if available. Usage of a python script
to pre-populate metadata files from the bio.tools api is described in
[Geting Started](../Getting_Started.md).

Metadata can also be obtained from tool manual and help pages,
[SciCrunch](https://scicrunch.org/), or other web resources.

## Common Metadata

The metadata fields for common metadata (common to whole tool suite)
files fields are described below.

### Required Fields

#### name

Name of the tool without subcommands. e.g. `grep`. Typically, extensions
are left out of the name. e.g. 'Picard', not 'picard.jar' although this
is not a requirement. Alternate names can optionally be stored in the
[alternateName](#alternatename) list.

#### identifier

A unique identifier for the tool. `capanno-utils` will generate this
automatically.

#### softwareVersion

The `softwareVersion` field has two subfields: `versionName` and
`indluededVersions`. `versionName` should follow the methods's version
convention as much as possible and can include variable portions,
e.g. v3.x.y. For some version conventions it may be appropriate for the
versionName to be a range, e.g. v301-v399. Exact supported versions
covered by the `versionName` can be specified in `includedVersions`.

Example:

``` yaml
softwareVersion: 
  versionName: 2.x
  includedVersions:
    - 2.1.1
    - 2.1.2
    - 2.2.1
    - 2.2.2
```

#### featureList

A list that contains the names of the subtools for a tool. The names in
the list must correspond to `subtoolName` provided in the metadata file
for the subtool, if one has been created, and should be the string used
to identify the subtool when running it from the command line, excluding
any dashes. If the tool can be called without a subcommand, the name of
the main tool is `__main__` and should be included in the list.

e.g. for pip:

``` yaml
featureList:
  - install
  - download
  - uninstall
  - freeze
  - list
  - show
  - check
  - config
  - search
  - wheel
  - hash
  - completion
  - help  # This one is not really necessary
```

#### metadataStatus

The status of the metadata file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

### Recommended Fields

#### description

Description of the method. Should be taken directly from the tool
documentation if it exists.

ex:

``` yaml
description: |
   grep  searches  for  PATTERN  in  each  FILE.  A FILE of “-” stands for standard input.  If no FILE is given, 
   recursive searches examine the working directory, and nonrecursive searches read standard input.  By default, 
   grep prints the matching lines.

   In addition, the variant programs egrep, fgrep and rgrep are the same as grep -E, grep -F, and grep -r, respectively.  
   These variants are deprecated, but are provided for backward compatibility.
```

#### codeRepository

Code repository information.

e.g. for `cwltool`

``` yaml
codeRepository:
  name: GitHub
  URL: https://github.com/common-workflow-language/cwltool
```

#### license

Software license for the tool defined by [SPDX
identifier](https://spdx.org/licenses/). e.g. `license: Apache-2.0`

#### WebSite

A list of websites associated with the tool, excluding the code
repository which should be provided in `codeRepository`. e.g.:

``` yaml
Website:
  - name: Samtools  # The name of the website.
    description: Samtools homepage.
    URL: http://www.htslib.org/
  - name: Wikipeda
    description: Wikipedia entry for SAMtools
    URL: https://en.wikipedia.org/wiki/SAMtools
```

#### contactPoint

A list that contains information about the software maintainer(s). Might
be convenient to add an anchor (with '&' symbol) to this field if the
maintainer(s) is also the also the tool creator so they can also be
referenced in the optional `creator` field. `identifier` fields are not
required, but if included should be a unique identifier such as orcid.

``` yaml
contactPoint: 
  - &Jane
    name: Jane Schmoe 
    email: jane@theschmoes.com
    identifier: https://orcid.org/0000-0001-6022-9825
  - name: Joe Schmoe
    email: joe@theschmoes.com
```

#### publication

A list of publications related to the tool. The first entry should be
the main reference to cite the tool.

e.g.

``` yaml
publication:
  - headline: 'BEDTools: a flexible suite of utilities for comparing genomic features.' # Title goes here.
    identifier: 10.1093/bioinformatics/btq033  # DOI goes here
  - headline: 'BEDTools: The Swiss-Army Tool for Genome Feature Analysis.'
    identifier: 10.1002/0471250953.bi1112s47
```

#### keywords

List of keywords that can be used as filters or tags to categorize the
method by topic and operation specified with an
[EDAM](http://bioportal.bioontology.org/ontologies/EDAM?p=classes) or
other ontology identifier. EDAM keywords are preferred. If you wish to
provide a keyword that is not in an ontology, it may be specified with
the keys `name` and `category`. The value of `category` must be either
'topic' or 'operation' for tools. These categories have the meanings
defined by [EDAM](http://edamontology.org/page).

``` yaml
keywords:
  - http://edamontology.org/operation_3182  # Genome alignment
  - http://edamontology.org/topic_0085  # Functional Genomics
  - name: ENCODE
    category: topic
```

### Optional Fields

#### alternateName

List of alternate names for the tool. This is convenient place to put
any other names for the tool that someone might use to search for it.

e.g. for `gzip`

``` yaml
alternateName: 
  - gunzip
  - zcat
```

#### creator

List of the tool's creator(s). Might be redundant if this is captured in
`contactPoint` and/or `publication`. Can alias `contactPoint` if this is
the case.

``` yaml
creator:
  - *Jane # alias to &Jane in the contactPoint field.
  - name: Bob Bobbins
    email: bob@bobsbobbins.eu
    identifier: 
```

#### programmingLanguage

List of programming languages that the tool is written in
`programmingLanguage: [Python2, C]`

#### datePublished

The release date of the software in international standard date notation
(YYYY-MM-DD). ex: `datePublished: 2003-04-14`

## Subtool metadata files

Metadata that is specific to subtools. If values are provided in the
subtool metadata and the common metadata, the value in subtool metadata
is used.

### Required Fields

#### name

The name must correspond to the name of the subtool as specified in the
[featureList](#featurelist) field of the common metadata file.

#### identifier

Unique identifier. This is generated automatically by `capanno-utils`.
The identifier format differs slightly for tools, subtools, scripts,
workflows, and use cases of these.

#### metadataStatus

The status of the subtool metadata file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### cwlStatus

The status of the subtool cwl file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### nextflowStatus

The status of the subtool nextflow file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### snakemakeStatus

The status of the subtool snakemake file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### wdlStatus

The status of the subtool wdl file.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

### Recommended Fields

#### description

Description that is specific to the subtool. Should be taken from the
documentation if available.

#### keywords

List of keywords that can be used as filters or tags to categorize the
method by topic and operation specified with an
[EDAM](http://bioportal.bioontology.org/ontologies/EDAM?p=classes) or
other ontology identifier. EDAM keywords are preferred. If you wish to
provide a keyword that is not in an ontology, it may be specified with
the keys `name` and `category`. The value of `category` must be either
'topic' or 'operation' for tools. These categories have the meanings
defined by [EDAM](http://edamontology.org/page).

``` yaml
keywords:
  - http://edamontology.org/operation_3182  # Genome alignment
  - http://edamontology.org/topic_0085  # Functional Genomics
  - name: ENCODE
    category: topic
```

### Optional Fields

#### alternateName

List of alternate names for the tool. This is convenient place to put
any other names for the tool that someone might use to search for it.

e.g. for `gzip`

``` yaml
alternateName: 
  - gunzip
  - zcat
```

# CWL Files

## Shared requirements

If multiple subtools or scripts have the same requirements
(SchemaDefRequirement, DockerRequirement, etc.) the requirement should
be moved to its own file, placed in the tool or script's common
directory, and be imported (using \$import) into each subtool/script CWL
file rather than being specified separately in each subtool/script CWL
file.

## Parameter completeness

The CWL documents in this repository are intended to describe use cases
for tools, subtools, and scripts. Documents should be in progress
towards including all the parameters for a method. The --help parameter
or other similar options are an exception. Authors are encouraged to
contribute a tools/subtools with incomplete parameter descriptions that
can be built upon later by themselves or others.

## File names

Files must be named with the format `{toolName}-{subtoolName}.cwl` or
`{sciptName}.cwl` The subtoolName must be excluded if the tool is called
without a subcommand.

## Defining name and versions

The name and version of the methods must be specified in the tool's
metadata.

## Field specific instructions

### label

May be used, but should be limited to a single line.

### doc

May be used to make notes about the CWL file. Documentation for the tool
itself should be placed in a separate metadata file.

### cwlVersion

Must be `v1.0`

### class

Must be `CommandLineTool`

### baseCommand and arguments

Command parameters that are ALWAYS present at the beginning of the
command, with immutable order must be placed in `baseCommand`. Other
arguments that do not correspond to input parameters must be placed in
`arguments`. [Picard tools](https://broadinstitute.github.io/picard/)
provides an example where arguments are needed because java arguments
may or may not be present: `java jvm-args -jar picard.jar picard-args`.
In this case `java` would be the baseCommand `jvm-args`, `-jar`,
`picard.jar` and `picard-args` would all be specified in `arguments`.
NOTE: `jvm-args` would also be described in `inputs` but would not have
an `inputBinding` field, since the command line binding instructions
would be described `arguments`.

### requirements and hints

Requirements needed to validate a CWL file must be placed in
`requirements`. All other 'requirements' must be placed in `hints`.

#### requirements that go in `requirements`

-   SchemaDefRequirement
-   InlineJavascriptRequirement

#### requirements that go in `hints`

-   DockerRequirement
-   SoftwareRequirement
-   InitialWorkDirectoryRequirement
-   EnvironmentVarRequirement
-   ResourceRequirement

### Tool inputs and outputs

#### types

-   Do not use `type: string` parameters for names of input or reference
    files/directories; use `type: File` or `type: Directory` as
    appropriate.

-   Use `type: enum` instead of `type: string` for elements with a fixed
    list of valid values.

-   Custom types should be defined with one external YAML per type
    definition for re-use.

-   Array and optional inputs may be specified with the shorthand
    syntaxes, `[]`, and `?` respectively, where allowed.

#### names

-   All input parameters that have an `inputBinding.prefix` field must
    be named the same as the prefix, minus any prepended dashes,
    maintaining any capitalization, internal dashes, and underscores.
    Abbreviated versions of the prefix should be avoided in favor of the
    long version. e.g. `-all`, not `-a`.

-   All `input` and `output` identifiers should be named as they are
    described in the tool's documentation. If documented names are
    unavailable, use informative names such as `unaligned_sequences`
    instead of generic names such as `input4`.

#### Other inputs/outputs fields

-   Each input and output parameter should contain a `doc` field. If
    available, the `doc` field should contain the parameter description
    as it appears in the documentation for the tool (man page, output of
    help, etc.).

-   `format` should be specified for all input and output `File`s.
    Bioinformatics tools should use format identifiers from
    [EDAM](http://edamontology.org/page).

-   Mark all input and output `File`s that are read from or written to
    in a streaming compatible way (only once, no random-access), as
    `streamable: true`.

-   Only use JavaScript expressions when necessary or when it greatly
    reduces complexity.

#### output glob

Be as specific as possible when describing glob patterns for outputs.
Authors should avoid wildcards if possible.

### namespaces and schemas

Include the EDAM namespace and schema.

``` yaml
$namespaces:
  edam: http://edamontology.org/
$schemas:
 - http://edamontology.org/EDAM_1.21.owl
```
