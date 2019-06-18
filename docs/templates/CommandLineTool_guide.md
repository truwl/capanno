# Command Line Tool Guide
## Introduction
This document contains xD Bio's best practices for the describing command line tools with the
[Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language) (CWL). 
Workflows, along with tool and workflow input files (job files) will be addressed in separate documents. 
Several practices here were adapted from 
[Common Workflow Language User Guide: Recommended Practices](http://www.commonwl.org/user_guide/rec-practices/)

### Definitions
The words "MUST", "REQUIRED",  "SHOULD" and "MAY" have the meanings described in 
[Key words for use in RFCs to Indicate Requirement Levels](https://www.ietf.org/rfc/rfc2119.txt)


### Scope of CommandLineTool CWL files

Tools that contain a limited number of arguments that can be provided together should be described by a single cwl file. 
Tools that contain multiple subtools or modes that include mutually exclusive arguments should be divided into 
separate CWL files. For example md5sum should be divided into `md5sum.cwl` and `md5sum-check.cwl` since `md5sum --check [File]` 
has arguments that are not relevant when running `md5sum`  without the --check option.

### CWL tool status
CWL tools may be contributed at any point of their development. 
The state of the CWL tool or subtool must be documented by using 
[semantic versioning](https://semver.org/spec/v2.0.0.html) and specified in the [tool's metadata](CommandLineTool_metadata_guide.md#version1). 


### Shared requirements
If multiple subtools have same requirements (SchemaDefRequirement, DockerRequirement, etc.) the requirement should be
moved to its own file, placed in the tool's common directory, and be imported (using $import) into each subtool CWL file 
rather than being specified separately in each subtool CWL file. 

### Parameter completeness
The cwl documents in this repository are intended to completely describe
use cases for tools and subtools. Documents should be in progress towards including all the parameters for a 
tool/subtool. The --help parameter or other similar options are an exception to this, although 
--help could be described as its own subtool. Authors are encouraged to contribute a tools/subtools with
incomplete parameter descriptions that can be built upon later by themselves or others.

### File names
Files must be named with the format `{toolName}-{subtoolName}.cwl` The subtoolName must be 
excluded if the tool is not divided into subtools. 
More on file names can be found in this repository's [README](../../README.md)

### Defining name and versions
The name and version of the tool must be specified in the [tool's metadata](CommandLineTool_metadata_guide.md)

## Field specific instructions
### label
May be used, but should be limited to a single line.
### doc
May be used to make notes about the CWL file.  Documentation for the tool itself should be placed in 
a separate metadata file as described in the [metadata guide](CommandLineTool_metadata_guide.md).

### cwlVersion
Must be v1.0

### baseCommand and arguments
Command parameters that are ALWAYS present at the beginning of the command, with immutable order
must be placed in `baseCommand`. Other arguments that do not correspond to input parameters must be
placed in `arguments`. [Picard tools](https://broadinstitute.github.io/picard/) provides an example where arguments are 
needed because java  arguments may or may not be present: `java jvm-args -jar picard.jar picard-args`. In this case `java` would
be the baseCommand `jvm-args`, `-jar`, `picard.jar` and  `picard-args` would all be specified in `arguments`. NOTE: 
`jvm-args` would also be described in `inputs` but would not have an `inputBinding` field, since the command line
binding instructions would be described `arguments`. 

### requirements and hints
Requirements needed to validate a tool must be placed in `requirements`. All other 'requirements'
must be placed in `hints`.
#### requirements that go in `requirements`
- SchemaDefRequirement
- InlineJavascriptRequirement

#### requirements that go in `hints`
- DockerRequirement
- SoftwareRequirement
- InitialWorkDirectoryRequirement
- EnvironmentVarRequirement
- ResourceRequirement

### Tool inputs and outputs

#### types

- Do not use `type: string` parameters for names of input or reference
files/directories; use `type: File` or `type: Directory` as appropriate.

- Use `type: enum` instead of `type: string` for elements with a fixed
list of valid values.

- Custom types should be defined with one external YAML per type
definition for re-use.

- Array and optional inputs may be specified with the shorthand syntaxes, `[]`, and `?` respectively, where allowed.

#### names
- All input parameters that have an `inputBinding.prefix` field must be 
named the same as the prefix, minus any prepended dashes, maintaining any 
capitalization, internal dashes, and underscores. Abbreviated versions of 
the prefix should be avoided in favor of the long version. e.g. `-all`, 
not `-a`.

- All `input` and `output` identifiers should be named as they are described in the tool's documentation. 
If documented names are unavailable, use informative names such as `unaligned_sequences` instead of generic names
 such as `input4`.

#### Other inputs/outputs fields 
- Each input and output parameter should contain a `doc` field. If available, 
the `doc` field should contain the parameter description as it appears in 
the documentation for the tool (man page, output of help, etc.).

- `format` should be specified for all input and output `File`s.
Bioinformatics tools should use format identifiers from [EDAM][edam-example].

- Mark all input and output `File`s that are read from or written to in a
streaming compatible way (only once, no random-access), as `streamable: true`.

- Only use JavaScript expressions when necessary or when it greatly reduces complexity. 

#### output glob
Be as specific as possible when describing glob patterns for outputs. Authors should
avoid wildcards if possible.

### namespaces and schemas
Include the EDAM namespace and schema.
~~~yaml
$namespaces:
  edam: http://edamontology.org/
$schemas:
 - http://edamontology.org/EDAM_1.21.owl
~~~

### Metadata
Metadata must be included in a separate metadata file. Optionally, the metadata may be imported into the CWL  file.

