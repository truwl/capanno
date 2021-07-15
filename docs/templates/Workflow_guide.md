## File names

Files should keep the names given by the original workflow author and be specified int the workflow's metadata file.

## Field specific instructions

### label

May be used, but should be limited to a single line.

### doc

May be used to make notes about the CWL file.  Documentation for the workflow itself should be placed in 
a separate metadata file.

### cwlVersion

Must be `v1.0`

### class

Must be `Workflow`

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

### inputs

### outputs

### steps