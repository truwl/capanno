# Introduction

This document contains Truwl's best practices for the describing command line tools with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL) and tool metadata
files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

## Definitions

!include docs/components/definitions.md

## CWL and metadata types 

There are 2 types of CWL tool files and three types of tool metadata files:

-  **Complete tools:** 
Metadata and CWL files for a tool that is not separated into separate subtools. 

-  **Parent metadata** 
 Metadata (only) that is inherited by multiple subtools.
 
- **Subtools**
Metadata and CWL tool files that are specific to a subtool. 

Metadata files can be initialized (along with the proper directory structure) programmatically as described in [getting started.]

## Scope of CommandLineTool CWL files

Tools that contain a limited number of arguments that can be provided together should be described by a single CWL and metadata file. Tools that contain multiple subtools or modes that include mutually exclusive arguments should be divided into 
separate CWL files and metadata files as described above. For example md5sum should be divided into `md5sum.cwl` and `md5sum-check.cwl` since `md5sum --check [File]` 
has arguments that are not relevant when running `md5sum`  without the --check option.

## CWL tool status

!include docs/components/cwl_file_status.md

# CWL Files

!include docs/templates/CommandLineTool_guide.md

# Tool metadata files

This describes how to document metadata for command line tools described by the common workflow language (CWL) in this repository. 
Separate metadata files must be specified for each version of a tool (excluding patch versions). Metadata files are
specified in YAML and can be generated programmatically (see [getting started](../docs/Getting_Started.md). Most metadata keys are defined by [schema.org](https://schema.org/) vocabularies and can be adapted to be 
imported by CWL files.

## List order

Metadata fields that accept lists as values should be ordered by importance if applicable. 

## Metadata sources

We recommend populating metadata files from [bio.tools](https://bio.tools/) if available. Usage of a python script to pre-populate 
metadata files from the bio.tools api is described in [geting started].

Metadata can also be obtained from tool manual and help pages, [SciCrunch](https://scicrunch.org/), or other web resources.

## <a name="complete"><a/>Tool and subtool metadata files

Metadata file fields for a complete tool or subtool. For subtools, fields marked with a \* can be provided in the metadata file for the subtool or can be inherited from a [parent metadata]() file, and fields marked with \*\*  must be inherited from the parent metadata file and cannot be provided in the primary subtool metadata file.  If values are provided in the subtool metadata and the parent metadata, the value in subtool metadata is used.

### Required Fields

### applicationSuite (subtools only)

!include docs/components/applicationSuite_field.md

#### <a name="name1"><a/>name

The name of the tool. Typically, extensions are left out of the name. e.g. 'Picard', not 'picard.jar' although this is
not a requirement. Alternate names can optionally be stored in the [alternateName](#alternatename) list. 
ex: `name: Picard`

For subtools, the name must correspond to the name of the subtool as specified in the [featureList](#featurelist) field of the primary tool metadata file.

ex: `search`

#### <a name="softwareVersion"></a> \*softwareVersion

!include docs/components/softwareVersion_field.md


#### <a name="version"></a> version

!include docs/components/version.md

#### identifier

!include docs/components/identifier_field.md


### Recommended Fields

#### \*description

!include docs/components/description_field.md

#### \*\*codeRepository

!include docs/components/codeRepository_field.md

#### \*\*license

!include docs/components/license_field.md

#### \*\*contactPoint

!include docs/components/contactPoint_field.md

#### \*\*publication

!include docs/components/publication_field.md

#### \*keywords

!include docs/components/keywords_field.md


### Optional Fields

#### \*alternateName

!include docs/components/alternateName_field.md

#### \*\*creator

!include docs/components/creator_field.md

#### \*\*programmingLanguage

!include docs/components/programmingLanguage_field.md

#### \*\*datePublished

!include docs/components/datePublished_field.md


## <a name="parent"><a/>Parent (common) Metadata

Metadata file fields for metadata that is inherited by multiple subtools. Fields that have already been described above are listed without a description.

### Required Fields

#### name

Name of the main tool. ex: `pip`.

#### softwareVersion

#### featureList

!include docs/components/featureList_field.md


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

