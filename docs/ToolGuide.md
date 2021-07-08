Introduction
============

This document contains Truwl's best practices for the describing
command line tools with metadata files and the [Common Workflow
Language](https://github.com/common-workflow-language/common-workflow-language)
(CWL). Several practices here were adapted from
[CommonsWorkflow Language User Guide: Recommended
Practices](http://www.commonwl.org/user_guide/rec-practices/). Best practices for other workflow language files (WDL, Snakemake, Nextflow) will be added as developed.

Definitions
-----------

The words "MUST", "REQUIRED", "SHOULD" and "MAY" have the meanings
described in [Key words for use in RFCs to Indicate Requirement
Levels](https://www.ietf.org/rfc/rfc2119.txt)

Division of tools into subtools
----------------------

There two types of tool metadata files:

-   **common-metadata** Metadata that is common to all subtools of .

-   **Subtools** Metadata and CWL tool files that are specific to a
    subtool.

Metadata files can be initialized (along with the proper directory
structure) programmatically as described in \[getting started.\]

Scope of CommandLineTool CWL files
----------------------------------

Tools that contain a limited number of arguments that can be provided
together should be described by a single CWL and metadata file. Tools
that contain multiple subtools or modes that include mutually exclusive
arguments should be divided into separate CWL files and metadata files
as described above. For example md5sum should be divided into
`md5sum` and `md5sum-check` since `md5sum --check [File]` has
arguments that are not relevant when running `md5sum` without the
--check option.

CWL tool status
---------------

CWL files may be contributed at any point of their development. The
state of the CWL file must be documented by using [semantic
versioning](https://semver.org/spec/v2.0.0.html) and specified in the
tool's metadata `version` field.

CWL Files
=========

Shared requirements
-------------------

If multiple subtools or scripts have the same requirements
(SchemaDefRequirement, DockerRequirement, etc.) the requirement should
be moved to its own file, placed in the tool or script's common
directory, and be imported (using \$import) into each subtool/script CWL
file rather than being specified separately in each subtool/script CWL
file.

Parameter completeness
----------------------

The cwl documents in this repository are intended to completely describe
use cases for tools, subtools, and scripts. Documents should be in
progress towards including all the parameters for a method. The --help
parameter or other similar options are an exception to this, although
--help could be described as its own subtool. Authors are encouraged to
contribute a tools/subtools with incomplete parameter descriptions that
can be built upon later by themselves or others.

File names
----------

Files must be named with the format `{toolName}-{subtoolName}.cwl` or
`{sciptName}.cwl` The subtoolName must be excluded if the tool is not
divided into subtools. More on file names can be found in this
repository's [README](../../README.md)

Defining name and versions
--------------------------

The name and version of the methods must be specified in the tool's
metadata.

Field specific instructions
---------------------------

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
    long version. e.g. `-all`, not `-a`.

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
    \[EDAM\]\[edam-example\].

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

``` {.yaml}
$namespaces:
  edam: http://edamontology.org/
$schemas:
 - http://edamontology.org/EDAM_1.21.owl
```

Tool metadata files
===================

This describes how to document metadata for command line tools described
by the common workflow language (CWL) in this repository. Separate
metadata files must be specified for each version of a tool (excluding
patch versions). Metadata files are specified in YAML and can be
generated programmatically (see [getting
started](../docs/Getting_Started.md). Most metadata keys are defined by
[schema.org](https://schema.org/) vocabularies and can be adapted to be
imported by CWL files.

List order
----------

Metadata fields that accept lists as values should be ordered by
importance if applicable.

Metadata sources
----------------

We recommend populating metadata files from
[bio.tools](https://bio.tools/) if available. Usage of a python script
to pre-populate metadata files from the bio.tools api is described in
\[geting started\].

Metadata can also be obtained from tool manual and help pages,
[SciCrunch](https://scicrunch.org/), or other web resources.

<a name="complete"><a/>Tool and subtool metadata files
------------------------------------------------------

Metadata file fields for a complete tool or subtool. For subtools,
fields marked with a \* can be provided in the metadata file for the
subtool or can be inherited from a [parent metadata]() file, and fields
marked with \*\* must be inherited from the parent metadata file and
cannot be provided in the primary subtool metadata file. If values are
provided in the subtool metadata and the parent metadata, the value in
subtool metadata is used.

### Required Fields

### applicationSuite (subtools only)

Information to identify the main tool. The subtool will inherit metadata
from the parent metadata. The `name` and `SoftwareVersion` fields must
match the corresponding fields in the \[main tool metadata file\].(\#)

``` {.yaml}
applicationSuite:
  name: pip
  softwareVersion: v19.0
  identifier:  # truwl identifier of main tool metadata, if it exists.
```

#### <a name="name1"><a/>name

The name of the tool. Typically, extensions are left out of the name.
e.g. 'Picard', not 'picard.jar' although this is not a requirement.
Alternate names can optionally be stored in the
[alternateName](#alternatename) list. ex: `name: Picard`

For subtools, the name must correspond to the name of the subtool as
specified in the [featureList](#featurelist) field of the primary tool
metadata file.

ex: `search`

#### <a name="softwareVersion"></a> \*softwareVersion

A string that specifies the version of the software that the CWL file is
valid for. Versions must follow the version conventions used by the
method author. ex: `softwareVersion: v3.2`

#### <a name="version"></a> version

Specifies the version of the CWL file and is used to determine its
development state. Must follow [semantic
versioning](https://semver.org/spec/v2.0.0.html) conventions. A 1.0 or
greater version of a CWL document must contain valid CWL syntax, follow
the required [best practices](CommandLineTool_guide.md), and describe
enough input parameters to be useful to others. Breaking changes that
constitute a major revision include changing the name/id of a command
input parameter that would make any job files not run properly. If not
specified, will be initialized to 0.1.0 ex: `version: 0.2.1`

#### identifier

Unique identifier for truwl.com. If not provided, it will be initialized
for you.

### Recommended Fields

#### \*description

Description of the method. Should be taken directly from the tool
documentation if it exists.

ex:

``` {.yaml}
description: |
   grep  searches  for  PATTERN  in  each  FILE.  A FILE of “-” stands for standard input.  If no FILE is given, 
   recursive searches examine the working directory, and nonrecursive searches read standard input.  By default, 
   grep prints the matching lines.

   In addition, the variant programs egrep, fgrep and rgrep are the same as grep -E, grep -F, and grep -r, respectively.  
   These variants are deprecated, but are provided for backward compatibility.
```

#### \*\*codeRepository

Code repository information.

ex: cwltool

``` {.yaml}
codeRepository:
  name: GitHub
  URL: https://github.com/common-workflow-language/cwltool
```

#### \*\*license

Software license for the tool defined by [SPDX
identifier](https://spdx.org/licenses/). e.g. `license: Apache-2.0`

#### \*\*contactPoint

A list that contains information about the software maintainer(s). Might
be convenient to add an anchor (with '&' symbol) to this field if the
maintainer(s) is also the also the tool creator so they can also be
referenced in the optional `creator` field. `identifier` fields are not
required, but if included should be a unique identifier such as orcid.

``` {.yaml}
contactPoint: 
  - &Jane
    name: Jane Schmoe 
    email: jane@theschmoes.com
    identifier: https://orcid.org/0000-0001-6022-9825
  - name: Joe Schmoe
    email: joe@theschmoes.com
```

#### \*\*publication

A list that describes publications related to the tool. The first entry
should be the main reference to cite the tool.

ex:

``` {.yaml}
publication:
  - headline: 'BEDTools: a flexible suite of utilities for comparing genomic features.' # Title goes here.
    identifier: 10.1093/bioinformatics/btq033  # DOI goes here
  - headline: 'BEDTools: The Swiss-Army Tool for Genome Feature Analysis.'
    identifier: 10.1002/0471250953.bi1112s47
```

#### \*keywords

List of tags to categorize the method by topic and operation specified
with an
[edam](http://bioportal.bioontology.org/ontologies/EDAM?p=classes) or
other ontology identifier. EDAM keywords are preferred. If you wish to
provide a keyword that is not in an ontology, it may be specified with
the keys `name` and `category`. The value of `category` must be either
'topic' or 'operation'. These categories have the meanings defined by
[EDAM](http://edamontology.org/page).

ex.

``` {.yaml}
keywords:
  - http://edamontology.org/operation_3182  # Genome alignment
  - http://edamontology.org/topic_0085  # Functional Genomics
  - name: ENCODE
    category: topic
```

### Optional Fields

#### \*alternateName

List of alternate names for the tool. This is convenient place to put
any other names for the tool that someone might use to search for it.
ex: `alternateName: []`

#### \*\*creator

List of the tool's creator(s). Might be redundant if this is captured in
`contactPoint` and/or `publication`. Can alias `contactPoint` if this is
the case.

``` {.yaml}
creator:
  - *Jane # alias to &Jane in the contactPoint field.
  - name: Bob Bobbins
    email: bob@bobsbobbins.eu
    identifier: 
```

#### \*\*programmingLanguage

List of programming languages that the tool is written in
`programmingLanguage: [Python2, C]`

#### \*\*datePublished

The release date of the particular version of software in international
standard date notation (YYYY-MM-DD). ex: `datePublished: 2003-04-14`

<a name="parent"><a/>Parent (common) Metadata
---------------------------------------------

Metadata file fields for metadata that is inherited by multiple
subtools. Fields that have already been described above are listed
without a description.

### Required Fields

#### name

Name of the main tool. ex: `pip`.

#### softwareVersion

#### featureList

Required for tools that are divided into subtools. The names in the list
must correspond to `subtoolName` provided in the metadata file for the
subtool, if one has been created, and should be the string used to
identify the subtool when running it from the command line, excluding
any dashes. `featureList` must include all the subtools of the tool even
when a metadata file or CWL file has not been created for the subtool.
When the `featureList` field is populated the metadata file will not
correspond to a single CWL file but will be used as a common metadata
file for multiple subtools to inherit from and must be placed in the
tool's `common/` directory.

ex. For pip:

``` {.yaml}
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

### Recommended Fields

#### description

#### codeRepository

#### codeRepository

#### license

#### WebSite

#### contactPoint

#### publication

#### keywords

### Optional Fields

#### alternateName

#### creator

#### programmingLanguage

#### datePublished
