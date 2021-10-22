## Introduction

This document contains Truwl's practices for describing worklfow metadata and workflow language files.

## Definitions

!include docs/components/definitions.md


## Workflow metadata files
Workflow metadata is specified in YAML files and can be initiated with the proper format and directory structure using `capanno-utils` as described in [Getting Started](./GettingStarted.md). Workflow descriptions differ from the how [tools](ToolGuide.md) are described in that only one main workflow language file is expected per workflow, and it will be in a single workflow language. A tool is broken up into various subtools and can be described separately by multiple workflow languages. Neither of these generally apply to workflows. It would be difficult to express the *same* workflow in multiple workflow languages and workflows are not typically written run in many diverse ways.

### Required Fields

#### name

The name of the workflow

#### softwareVersion

!include docs/components/softwareVersion_field.md

#### identifier

!include docs/components/identifier_field.md

#### metadataStatus

The status of the workflow's metadata.

!include docs/components/status_field.md

#### workflowStatus

The status of the workflow.

!include docs/components/status_field.md

#### workflowLanguage

The language that the workflow is written in. Must be one of cwl, wdl, nextflow, or snakemake.

#### workflowFile

The name of the main workflow file. e.g. `atac.wdl`

#### repoName

The name of the git repository if public.

#### gitTag

The git tag for the repo if available; usually a version.

#### inputsTemplate

The name of an inputs file template: a **blank** inputs.json for WDL, job.yaml for CWL, etc.

#### callMap

!include docs/components/callMap_field.md

### Recommended Fields

#### description

!include docs/components/description_field.md

#### codeRepository

!include docs/components/codeRepository_field.md

#### license

!include docs/components/license_field.md

#### contactPoint

!include docs/components/contactPoint_field.md

#### publication

!include docs/components/publication_field.md

#### keywords

!include docs/components/keywords_field.md

### Optional Fields

#### alternateName

!include docs/components/alternateName_field.md

#### creator

!include docs/components/creator_field.md

#### datePublished

!include docs/components/datePublished_field.md



## CWL Files

!include docs/templates/Workflow_guide.md

## WDL Files

No guidance at this point.

## Snakemake Files

No guidance at this point.

## Nextflow Files

No guidance at this point.