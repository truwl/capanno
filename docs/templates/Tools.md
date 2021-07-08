# Introduction

This document contains Truwl's best practices for the describing command line tools with metadata files and
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL)
files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/). Best practices for describing workflow language files from other workflow languages will be added as developed.

## Definitions

!include docs/components/definitions.md

## Dividing tools into subtools
All tools in `capanno` are divided into subtools. Subtools are the subcommands or modes that are specified when the tool is called e.g. `<tool_name> <subtool_name> <arguments>`. For tools that can be called without a subcommand, the subtool name is `__main__` and the files associated with the main tool are placed in a directory with the same named after the tool name, e.g. of for cat `tools/cat/8.x/cat`. 

### Metadata files
There are 2 types of metadata files for tools:

-  **Common metadata**: Metadata that is common to all subtools. The main source of metadata for all subtools is typically from common metadata.
 
- **Subtool metadata**:
Metadata in subtool metadata files can add metadata specific to the subtool or override metadata from the common metadata.

Metadata files can be initialized (along with the proper directory structure) programmatically as described in [Getting Started](../Getting_Started.md).

### Scope of tool workflow language files

Workflow language files must be divided into the same subtool structure as metadata, i.e., there must be only one workflow language file per subtool, per workflow language, i.e. a single subtool may have one of each CWL, WDL, Snakemake, and Nextflow files. Subtools typically contain mutually exclusive arguments that should be divided into 
separate workflow language files. For example md5sum CWL file would be divided into `md5sum/<versionName>/md5sum/md5sum.cwl` and `md5sum/<versionName>/md5sum_check/md5sum-check.cwl` since `md5sum --check [File]` 
has arguments that are not relevant when running `md5sum`  without the --check subcommand.


# Tool metadata files

This describes how to document metadata for command line tools described in this repository. Metadata files are
specified in YAML and can be generated programmatically (see [getting started](../../docs/Getting_Started.md). Most metadata keys are defined by [schema.org](https://schema.org/) vocabularies and can be adapted to be 
imported by CWL files.

## List order

Metadata fields that accept lists as values should be ordered by importance if applicable. 

## Metadata sources

We recommend populating metadata files from [bio.tools](https://bio.tools/) if available. Usage of a python script to pre-populate 
metadata files from the bio.tools api is described in [Geting Started](../Getting_Started.md).

Metadata can also be obtained from tool manual and help pages, [SciCrunch](https://scicrunch.org/), or other web resources.

## <a name="common"><a/>Common Metadata

The metadata fields for common metadata files fields are described below.

### Required Fields

#### name

Name of the tool without subcommands. e.g. `grep`.  Typically, extensions are left out of the name. e.g. 'Picard', not 'picard.jar' although this is
not a requirement. Alternate names can optionally be stored in the [alternateName](#alternatename) list.

#### identifier
A unique identifier for the tool. `capanno-utils` will generate this automatically.

#### softwareVersion

!include docs/components/softwareVersion_field.md

#### <a name="featureList"><a/> featureList

!include docs/components/featureList_field.md

#### metadataStatus
The status of the metadata file.
!include docs/components/status_field.md

### Recommended Fields 

#### description

!include docs/components/description_field.md

#### codeRepository

!include docs/components/codeRepository_field.md

#### license

!include docs/components/license_field.md

#### WebSite

!include docs/components/license_field.md

#### contactPoint

!include docs/components/contactPoint_field.md

#### publication

!include docs/components/publication_field.md

#### keywords

!include docs/components/publication_field.md

### Optional Fields

#### <a name="alternate"><a/>alternateName

!include docs/components/alternateName_field.md

#### creator

#### programmingLanguage

#### datePublished



## <a name="subtools"><a/>Subtool metadata files

For subtools, fields marked with a * can be provided in the metadata file for the subtool or can be inherited from a [common metadata]() file, and fields marked with \*\*  must be inherited from the common metadata file and cannot be provided in the subtool metadata file.  If values are provided in the subtool metadata and the common metadata, the value in subtool metadata is used.

### Required Fields

#### <a name="name1"><a/>name

The name must correspond to the name of the subtool as specified in the [featureList](#featurelist) field of the common metadata file.


#### <a name="version"></a> version
No longer used for subtools.

#### identifier

!include docs/components/identifier_field.md

#### metadataStatus
The status of the subtool metadata file.
!include docs/components/status_field.md

#### cwlStatus
The status of the subtool cwl file.
!include docs/components/status_field.md

#### nextflowStatus
The status of the subtool nextflow file.
!include docs/components/status_field.md

#### snakemakeStatus
The status of the subtool snakemake file.
!include docs/components/status_field.md

#### wdlStatus
The status of the subtool wdl file.
!include docs/components/status_field.md

### Recommended Fields

#### description

!include docs/components/description_subtool.md

#### \*keywords

!include docs/components/keywords_field.md


### Optional Fields

#### \*alternateName

!include docs/components/alternateName_field.md

# CWL Files

!include docs/templates/CommandLineTool_guide.md