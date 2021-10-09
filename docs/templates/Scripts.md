
# Introduction

This document contains xD Bio's best practices for the describing scripts with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL) and script metadata files. Scripts are CWL CommandLineTools but are described differently than a regular tool in their metadata files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

## Definitions

!include docs/components/definitions.md

## CWL and metadata types 

There is a single type of CWL script file and two types of script metadata files:

-  **scripts:** 
Metadata and CWL files for a script,

-  **common script metadata** 
 Metadata (only) that is inherited by multiple scripts.
 

Metadata files can be initialized (along with the proper directory structure) programmatically as described in [getting started](../docs/Getting_Started.md)

## Status of script CWL files

!include docs/components/cwl_file_status.md

NOTE: script cwl files follow the same best practices as tool cwl files. The guide is repeated here for convenience.

# CWL Files

!include docs/templates/CommandLineTool_guide.md

# Script metadata files

## Script Metadata

Metadata for a script. Fields marked with a \* can be provided in the primary metadata file, or can be inherited from a common script metadata file.

### Required

#### <a name="name1"><a/>name

The name of the script. Alternate names can optionally be stored in the [alternateName](#alternatename) list. 
ex: `name: foo.py`

#### \*softwareVersion

!include docs/components/softwareVersion_field.md

#### identifier

!include docs/components/identifier_field.md

#### parentMetadata

Required for script metadata that inherits some of its data from a common file. 

!include docs/components/parentMetadata_field.md

### Recommended Fields

#### \*description

!include docs/components/description.md

#### \*codeRepository

!include docs/components/codeRepository_field.md

#### \*WebSite

!include docs/components/WebSite_field.md

#### \*license

!include docs/components/license_field.md

#### \*contactPoint

!include docs/components/contactPoint_field.md

#### \*publication

!include docs/components/publication_field.md

#### \*keywords

!include docs/components/keywords_field.md

#### \*parentScripts

!include docs/components/parentScripts_field.md

#### \*tools

!include docs/components/tools_field.md

### Optional Fields

#### alternateName

!include docs/components/alternateName_field.md

#### \*creator

!include docs/components/creator_field.md

#### \*programmingLanguage

!include docs/components/programmingLanguage_field.md

#### \*datePublished

!include docs/components/datePublished_field.md

## <a name="common"></a> Common Script Metadata

Metadata that can be inherited by multiple scripts.

Fields that can be placed in a common script metadata file and inherited by a regular script metadata file are marked with a \* in the Script Metadata descriptions.