# Introduction
This document contains xD Bio's best practices for the describing command line tools with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL) and tool metadata
files.
Several practices here were adapted from 
[CommonsWorkflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

## Definitions
The words "MUST", "REQUIRED",  "SHOULD" and "MAY" have the meanings described in 
[Key words for use in RFCs to Indicate Requirement Levels](https://www.ietf.org/rfc/rfc2119.txt)

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
CWL tools may be contributed at any point of their development. 
The state of the CWL tool or subtool must be documented by using 
[semantic versioning](https://semver.org/spec/v2.0.0.html) and specified in the tool's metadata `version` field.


# CWL Files

!include docs/templates/CommandLineTool_guide.md

# Tool metadata files

!include docs/templates/CommandLineTool_metadata_guide.md