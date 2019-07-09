## Introduction

This documents how to describe metadata for command line tools described by the common workflow language (CWL) in this repository. 
Separate metadata files must be specified for each version of a tool (excluding patch versions). Metadata files are
specified in YAML. Most metadata keys are defined by [schema.org](https://schema.org/) vocabularies and can be adapted to be 
imported by CWL files.

## List order

Metadata fields that accept lists as values should be ordered by importance if applicable. 

## Metadata sources

We recommend populating metadata files from [bio.tools](https://bio.tools/) if available. Usage of a python script to pre-populate 
metadata files from the bio.tools api is described in [geting started].

Metadata can also be obtained from tool manual and help pages, [SciCrunch](https://scicrunch.org/), or other web resources.

!include docs/components/version.md

## <a name="complete"><a/>Complete tool metadata files

Metadata file fields for a tool that is not separated into subtools.

### Required Fields

#### <a name="name1"><a/>name

The name of the tool. Typically, extensions are left out of the name. e.g. 'Picard', not 'picard.jar' although this is
not a requirement. Alternate names can optionally be stored in the [alternateName](#alternatename) list. 
ex: `name: Picard`

!include docs/components/softwareVersion_field.md


### Recommended Fields

!include docs/components/description_field.md

!include docs/components/codeRepository_field.md

!include docs/components/license_field.md

!include docs/components/contactPoint_field.md

!include docs/components/publication_field.md

!include docs/components/keywords_field.md


### Optional Fields

!include docs/components/alternateName_field.md

!include docs/components/creator_field.md

!include docs/components/programmingLanguage_field.md

!include docs/components/datePublished_field.md

!include docs/components/downloadURL_field.md


## <a name="parent"><a/>Parent (common) Metadata

Metadata file fields for metadata that is inherited by multiple subtools.

### Required Fields

#### name

Name of the main tool. ex: `pip`. Also see [name](#name1) above.

!include docs/components/softwareVersion_field.md

!include docs/components/featureList_field.md


### Recommended Fields

!include docs/components/description_field.md

!include docs/components/codeRepository_field.md

!include docs/components/license_field.md

!include docs/components/WebSite_field.md

!include docs/components/contactPoint_field.md

!include docs/components/publication_field.md

!include docs/components/keywords_field.md


### Optional Fields

!include docs/components/alternateName_field.md

!include docs/components/creator_field.md

!include docs/components/programmingLanguage_field.md

!include docs/components/datePublished_field.md

!include docs/components/downloadURL_field.md


## <a name="subtools"><a/>Subtool metadata

Metadata file fields for a CWL file that describes a subtool of a tool.

### Required Fields

!include docs/components/applicationSuite_field.md

#### name

This must correspond to the name of the subtool as specified in the [featureList](#featurelist) field of the primary tool metadata file.
ex: `search`

#### version

Same as [version](#version1) above. Each subtool CWL file
must have its own version.

### Recommended Fields

!include docs/components/description_field.md

Description of the subtool. Should contain subtool specific information.

ex. pip search

```yaml
description:  |
  Search PyPI for packages.
```

!include docs/components/keywords_field.md

Terms that are specfic to the subtool. 
Terms that are relevant to all subtools of a tool should be defined in the common metadata file.
See [keywords](#keywords) above.

### Optional Fields

!include docs/components/alternateName_field.md

List of alternate names for the subtool. 
ex: gzip -c . gzip is parent tool name. 'c' is subtool name. zcat is an alternateName.

```yaml
alternateName:
  - zcat
```
