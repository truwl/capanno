Introduction
============

This document contains xD Bio's best practices for the describing
command line tools with the [Common Workflow
Language](https://github.com/common-workflow-language/common-workflow-language)
(CWL) and tool metadata files. Several practices here were adapted from
[CommonsWorkflow Language User Guide: Recommended
Practices](http://www.commonwl.org/user_guide/rec-practices/)

Definitions
-----------

The words "MUST", "REQUIRED", "SHOULD" and "MAY" have the meanings
described in [Key words for use in RFCs to Indicate Requirement
Levels](https://www.ietf.org/rfc/rfc2119.txt)

CWL and metadata types
----------------------

There is one types of CWL workflow file and one type of workflow
metadata files.

CWL Workflow status
-------------------

CWL files may be contributed at any point of their development. The
state of the CWL file must be documented by using [semantic
versioning](https://semver.org/spec/v2.0.0.html) and specified in the
tool's metadata `version` field.

CWL Files
=========

Still lots to do here.

File names
----------

Files must be named with the format `{workflowName}.cwl`. Whitespace in
the workflow name should be replaced with underscores. More on file
names and directory structure can be found in this repository's
[README](../../README.md)

Defining name and versions
--------------------------

The name and version of the methods must be specified in the tool's
metadata.

Field specific instructions
---------------------------

### label

May be used, but should be limited to a single line.

### doc

May be used to make notes about the CWL file. Documentation for the
workflow itself should be placed in a separate metadata file.

### cwlVersion

Must be `v1.0`

### class

Must be `Workflow`

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

### inputs

### outputs

### steps

Workflow metadata files
=======================

### Required Fields

#### name

The name of the workflow

#### softwareVersion

A string that specifies the version of the software that the CWL file is
valid for. Versions must follow the version conventions used by the
method author. ex: `softwareVersion: v3.2`

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

#### identifier

Unique identifier for truwl.com. If not provided, it will be initialized
for you.

#### callMap

The `callMap` field relates the steps in a workflow to the tool, script,
or or workflow that is called. All steps must be in the callMap list.

ex:

``` {.yaml}
callMap:
  - id:  # id (name) of the workflow step
    identifier: # truwl identifier of the tool, script, or workflow that is called.
```

### Recommended Fields

#### description

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

#### codeRepository

Code repository information.

ex: cwltool

``` {.yaml}
codeRepository:
  name: GitHub
  URL: https://github.com/common-workflow-language/cwltool
```

#### license

Software license for the tool defined by [SPDX
identifier](https://spdx.org/licenses/). e.g. `license: Apache-2.0`

#### contactPoint

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

#### publication

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

#### keywords

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

#### alternateName

List of alternate names for the tool. This is convenient place to put
any other names for the tool that someone might use to search for it.
ex: `alternateName: []`

#### creator

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

#### programmingLanguage

List of programming languages that the tool is written in
`programmingLanguage: [Python2, C]`

#### datePublished

The release date of the particular version of software in international
standard date notation (YYYY-MM-DD). ex: `datePublished: 2003-04-14`
