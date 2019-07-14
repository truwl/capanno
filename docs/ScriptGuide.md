Introduction
============

This document contains xD Bio's best practices for the describing
scripts with the [Common Workflow
Language](https://github.com/common-workflow-language/common-workflow-language)
(CWL) and script metadata files. Scripts are CWL CommandLineTools but
are described differently than a regular tool in their metadata files.
Several practices here were adapted from [CommonsWorkflow Language User
Guide: Recommended
Practices](http://www.commonwl.org/user_guide/rec-practices/)

Definitions
-----------

The words "MUST", "REQUIRED", "SHOULD" and "MAY" have the meanings
described in [Key words for use in RFCs to Indicate Requirement
Levels](https://www.ietf.org/rfc/rfc2119.txt)

CWL and metadata types
----------------------

There is a single type of CWL script file and two types of script
metadata files:

-   **scripts:** Metadata and CWL files for a script,

-   **common script metadata** Metadata (only) that is inherited by
    multiple scripts.

Metadata files can be initialized (along with the proper directory
structure) programmatically as described in [getting
started](../docs/Getting_Started.md)

Status of script CWL files
--------------------------

CWL files may be contributed at any point of their development. The
state of the CWL file must be documented by using [semantic
versioning](https://semver.org/spec/v2.0.0.html) and specified in the
tool's metadata `version` field.

NOTE: script cwl files follow the same best practices as tool cwl files.
The guide is repeated here for convenience.

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

Script metadata files
=====================

Script Metadata
---------------

Metadata for a script. Fields marked with a \* can be provided in the
primary metadata file, or can be inherited from a common script metadata
file.

### Required

#### <a name="name1"><a/>name

The name of the script. Alternate names can optionally be stored in the
[alternateName](#alternatename) list. ex: `name: foo.py`

#### \*softwareVersion

A string that specifies the version of the software that the CWL file is
valid for. Versions must follow the version conventions used by the
method author. ex: `softwareVersion: v3.2`

#### identifier

Unique identifier for truwl.com. If not provided, it will be initialized
for you.

#### version

Specifies the version of the CWL file and is used to determine its
development state. Must follow [semantic
versioning](https://semver.org/spec/v2.0.0.html) conventions. A 1.0 or
greater version of a CWL document must contain valid CWL syntax, follow
the required [best practices](CommandLineTool_guide.md), and describe
enough input parameters to be useful to others. Breaking changes that
constitute a major revision include changing the name/id of a command
input parameter that would make any job files not run properly. If not
specified, will be initialized to 0.1.0 ex: `version: 0.2.1`

#### parentMetadata

Required for script metadata that inherits some of its data from a
common file.

Relative path (or list of paths) to metadata file(s) that contains
fields that the primary metadata should inherit from. Fields from
parentMetadata are only used if values are not provided in the primary
file and metaadata listed earlier takes precedence over metadata listed
later.

ex.

``` {.yaml}
parentMetadata: ../common/foo-metadata.yaml
```

or

``` {.yaml}
parentMetadata:
  - ../common/foo-metadata.yaml  # fields in foo-metadata.yaml take precedence to fields in bar-metadata.yaml
  - ../common/bar-metadata.yaml
```

### Recommended Fields

#### \*description

!include docs/components/description.md

#### \*codeRepository

Code repository information.

ex: cwltool

``` {.yaml}
codeRepository:
  name: GitHub
  URL: https://github.com/common-workflow-language/cwltool
```

#### \*WebSite

A list of websites associated with the tool, excluding the code
repository which should be provided in `codeRepository`.

ex:

``` {.yaml}
Website:
  - name: Samtools  # The name of the website.
    description: Samtools homepage.
    URL: http://www.htslib.org/
  - name: Wikipeda
    description: Wikipedia entry for SAMtools
    URL: https://en.wikipedia.org/wiki/SAMtools
```

#### \*license

Software license for the tool defined by [SPDX
identifier](https://spdx.org/licenses/). e.g. `license: Apache-2.0`

#### \*contactPoint

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

#### \*publication

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

#### \*parentScripts

List of scripts that the script uses (imports of other scripts). Used to
create navigable relationships.

ex.

``` {.yaml}
parentScripts:
  - name: 
    alternateName:  # useful for looking up scripts if searching 'name' does not provide hits.
    softwareVersion:  # Needed to find correct version of script.
    identifier: # truwl identifer if it exists. If this is provided, other fields not used.
```

#### \*tools

List of tools that the script calls. Used to create navigable
relationships.

ex.

``` {.yaml}
tools:
  - name:  #
    softwareVersion:
    identifier: # truwl identifer if it exists. If this is provided, other fields not used.
```

### Optional Fields

#### alternateName

List of alternate names for the tool. This is convenient place to put
any other names for the tool that someone might use to search for it.
ex: `alternateName: []`

#### \*creator

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

#### \*programmingLanguage

List of programming languages that the tool is written in
`programmingLanguage: [Python2, C]`

#### \*datePublished

The release date of the particular version of software in international
standard date notation (YYYY-MM-DD). ex: `datePublished: 2003-04-14`

<a name="common"></a> Common Script Metadata
--------------------------------------------

Metadata that can be inherited by multiple scripts.

Fields that can be placed in a common script metadata file and inherited
by a regular script metadata file are marked with a \* in the Script
Metadata descriptions.
