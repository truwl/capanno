# Introduction

This document contains xD Bio's best practices for the describing command line tools with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL) and tool metadata
files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

## Definitions

!include docs/components/definitions.md

## CWL and metadata types 

There is one types of CWL workflow file and one type of workflow metadata files.

## CWL Workflow status

!include docs/components/cwl_file_status.md

# CWL Files

Still lots to do here.

!include docs/templates/Workflow_guide.md

# Workflow metadata files

### Required Fields

#### name

The name of the workflow

#### softwareVersion

!include docs/components/softwareVersion_field.md

#### version

!include docs/components/version.md

#### identifier

!include docs/components/identifier_field.md

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

#### programmingLanguage

!include docs/components/programmingLanguage_field.md

#### datePublished

!include docs/components/datePublished_field.md
