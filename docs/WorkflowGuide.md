## Introduction

This document contains Truwl's practices for describing worklfow
metadata and workflow language files.

## Definitions

The words "MUST", "REQUIRED", "SHOULD" and "MAY" have the meanings
described in [Key words for use in RFCs to Indicate Requirement
Levels](https://www.ietf.org/rfc/rfc2119.txt)

## Workflow metadata files

Workflow metadata is specified in YAML files and can be initiated with
the proper format and directory structure using `capanno-utils` as
described in [Getting Started](./GettingStarted.md). Workflow
descriptions differ from the how [tools](ToolGuide.md) are described in
that only one main workflow language file is expected per workflow, and
it will be in a single workflow language. A tool is broken up into
various subtools and can be described separately by multiple workflow
languages. Neither of these generally apply to workflows. It would be
difficult to express the *same* workflow in multiple workflow languages
and workflows are not typically written run in many diverse ways.

### Required Fields

#### name

The name of the workflow

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

#### identifier

Unique identifier. This is generated automatically by `capanno-utils`.
The identifier format differs slightly for tools, subtools, scripts,
workflows, and use cases of these.

#### metadataStatus

The status of the workflow's metadata.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### workflowStatus

The status of the workflow.

Status can be set to `Incomplete`, `Draft`, or `Released`. Default =
`Incomplete`. Incomplete signifies that the file does not exist or has
been initalized but not populated. Draft signifies that the file is a
work in progress. Released means the file is (reasonably) ready to be
used.

#### workflowLanguage

The language that the workflow is written in. Must be one of cwl, wdl,
nextflow, or snakemake.

#### workflowFile

The name of the main workflow file. e.g. `atac.wdl`

#### repoName

The name of the git repository if public.

#### gitTag

The git tag for the repo if available; usually a version.

#### inputsTemplate

The name of an inputs file template: a **blank** inputs.json for WDL,
job.yaml for CWL, etc.

#### callMap

The `callMap` field relates the steps in a workflow to the tool, script,
or or workflow that is called. All steps must be in the callMap list.

ex:

``` yaml
callMap:
  - id:  # id (name) of the workflow step
    identifier: # truwl identifier of the tool, script, or workflow that is called.
```

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

#### datePublished

The release date of the software in international standard date notation
(YYYY-MM-DD). ex: `datePublished: 2003-04-14`

## CWL Files

## File names

Files should keep the names given by the original workflow author and be
specified int the workflow's metadata file.

## Field specific instructions

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

## WDL Files

No guidance at this point.

## Snakemake Files

No guidance at this point.

## Nextflow Files

No guidance at this point.
