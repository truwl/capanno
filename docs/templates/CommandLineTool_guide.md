## Shared requirements

If multiple subtools or scripts have the same requirements (SchemaDefRequirement, DockerRequirement, etc.) the requirement should be
moved to its own file, placed in the tool or script's common directory, and be imported (using $import) into each subtool/script CWL file 
rather than being specified separately in each subtool/script CWL file. 

## Parameter completeness

The cwl documents in this repository are intended to completely describe
use cases for tools, subtools, and scripts. Documents should be in progress towards including all the parameters for a 
method. The --help parameter or other similar options are an exception to this, although 
--help could be described as its own subtool. Authors are encouraged to contribute a tools/subtools with
incomplete parameter descriptions that can be built upon later by themselves or others.

## File names

Files must be named with the format `{toolName}-{subtoolName}.cwl` or `{sciptName}.cwl` The subtoolName must be 
excluded if the tool is not divided into subtools. 
More on file names can be found in this repository's [README](../../README.md)

## Defining name and versions

The name and version of the methods must be specified in the tool's metadata.

## Field specific instructions

### label

May be used, but should be limited to a single line.

### doc

May be used to make notes about the CWL file.  Documentation for the tool itself should be placed in 
a separate metadata file.

### cwlVersion

Must be `v1.0`

### class

Must be `CommandLineTool`

### baseCommand and arguments

Command parameters that are ALWAYS present at the beginning of the command, with immutable order
must be placed in `baseCommand`. Other arguments that do not correspond to input parameters must be
placed in `arguments`. [Picard tools](https://broadinstitute.github.io/picard/) provides an example where arguments are 
needed because java  arguments may or may not be present: `java jvm-args -jar picard.jar picard-args`. In this case `java` would
be the baseCommand `jvm-args`, `-jar`, `picard.jar` and  `picard-args` would all be specified in `arguments`. NOTE: 
`jvm-args` would also be described in `inputs` but would not have an `inputBinding` field, since the command line
binding instructions would be described `arguments`. 

### requirements and hints

Requirements needed to validate a CWL file must be placed in `requirements`. All other 'requirements'
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

```yaml
$namespaces:
  edam: http://edamontology.org/
$schemas:
 - http://edamontology.org/EDAM_1.21.owl
```
