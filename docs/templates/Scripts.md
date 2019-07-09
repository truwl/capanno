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
