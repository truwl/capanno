
# Introduction

This document contains xD Bio's best practices for the describing scripts with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL) and script metadata
files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

## Definitions

The words "MUST", "REQUIRED",  "SHOULD" and "MAY" have the meanings described in 
[Key words for use in RFCs to Indicate Requirement Levels](https://www.ietf.org/rfc/rfc2119.txt)

## CWL and metadata types 

There is a single type of CWL script file and two types of script metadata files:

-  **Scripts:** 
Metadata and CWL files for a script,

-  **common_script metadata** 
 Metadata (only) that is inherited by multiple scripts.
 

Metadata files can be initialized (along with the proper directory structure) programmatically as described in [getting started.]


## CWL tool status

CWL tools may be contributed at any point of their development. 
The state of the CWL tool or subtool must be documented by using 
[semantic versioning](https://semver.org/spec/v2.0.0.html) and specified in the tool's metadata `version` field.

NOTE: script cwl files follow the same best practices as tool cwl files. The guide is repeated here for conveinience.

# CWL Files

!include docs/templates/CommandLineTool_guide.md

# Script metadata files

## Script Metadata

Metadata for a script.

### Required

#### <a name="name1"><a/>name

The name of the script. Alternate names can optionally be stored in the [alternateName](#alternatename) list. 
ex: `name: foo.py`

!include docs/components/softwareVersion_field.md

#### parentMetadata

Required for script metadata that inherits some of its data from a common file. List of relative paths to [common script metadata files](#common) listed in order of specificity. Metadata listed earlier takes precedence over metadata listed later.
ex.

```yaml
parentMetadata:
  - ../common/myscripts_common1-metadata.yaml
  - ../common/myscripts_common2-metadata.yaml

```

### Recommended Fields

!include docs/components/description.md

!include docs/components/identifier_field.md

!include docs/components/version.md

!include docs/components/codeRepository_field.md

!include docs/components/WebSite_field.md

!include docs/components/license_field.md

!include docs/components/contactPoint_field.md

!include docs/components/publication_field.md

!include docs/components/keywords_field.md

#### parentScripts

List of scripts that the script uses. Used to create navigable relationships.

ex.

```yaml
parentScripts:
  - name: 
    alternateName:  # if looking 
    softwareVersion:
    identifier: # truwl identifer if it exists. If this is provided, other fields not used.

```

#### tools

List of tools that the script calls. Used to create navigable relationships.

ex.

```yaml
tools:
  - name:  #
    softwareVersion:
    identifier: # truwl identifer if it exists. If this is provided, other fields not used.

```

### Optional Fields

!include docs/components/alternateName_field.md

!include docs/components/creator_field.md

!include docs/components/programmingLanguage_field.md


!include docs/components/datePublished_field.md

## <a name="common"></a> Common Script Metadata

Metadata that can be inherited by multiple scripts.



