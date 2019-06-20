## Repository scope and purpose
This repository contains Common Workflow Language (CWL) files and associated metadata to describe and run 
command line tools and workflows with a an emphasis on tools used in bioinformatics. 


The goal is to create a set CWL tool wrappers that are thorough, well documented, consistent, 
openly curated, and widely used by the community.

identifiers.

Consistent tools, versioning, etc. Build workflows through tools.

instances.

Tools in this repository will be made available through [truwl.com](https://truwl.com).

## How to contribute
Contributions to CWL files, metadata files, best practices and other documentation files 
are highly encouraged and can be made through pull requests and issues. Cwl and metadata file contributions must follow
 the best practices described in  [best practices](docs/templates/CommandLineTool_guide.md). 

## Repository structure and directory/file names
This repository has three main content directories: cwl-tools, cwl-scripts, and cwl-workflows.
### cwl-tools
Description. What is a tool?

#### `cwl-tools/` subdirectories
The directories in `cwl-tools/` must be named with the name of the tool being described preserving capitalization, 
dashes, and underscores. 

### `cwl-tools/{toolName}/` subdirectories
The subdirectories for each tool must be named according to the version of the tool that is described
preserving the version conventions of the tool author. e.g. v3.2, 3.2, 3b, etc. Patch versions may be omitted.
Separate directories should be created for each version of the tool even if the contents of the 
directories are initially identical to allow the contents to diverge.

### `cwl-tools/{toolName}/{toolVersion}/` subdirectories 
Tools that contain a limited number of arguments that can be provided together should be described by a single cwl file
and metadata file that are contained in a directory named with the tool name. 
The CWL file must be named `{toolName}.cwl` 
and the metadata file must be named `{toolName}-metadata.yaml`. 
e.g. `cwl-tools/cat/8.25/cat/cat.cwl` and `cwl-tools/cat/8.25/cat/cat-metadata.yaml` 

Tools that contain multiple subtools or modes should be divided into rational components (subtools) and be placed in their 
own directory named `{toolName}-{subtoolName}`.  
e.g. `cwl-tools/STAR/2.7/STAR-genomeGenerate/STAR-genomeGenerate.cwl` and `cwl-tools/STAR/2.7/STAR-genomeGenerate/STAR-genomeGenerate-metadata.yaml`
If dividing subtools into more finite components is deemed necessary, more directory nesting must not be added and the further breakdown
must be specified in the filename. 
e.g. `cwl-tools/{tool_name}/{tool_version}/{tool_name}-{subtool_name}/{tool_name}-{subtool_name}-{subtool_component_name}.cwl`
The metadata file for a subtool should be specific to the subtool. 
Metadata pertinent to all subtools should be placed in a parent metadata file in the  the `common/` directory described below. 
If particular metadata fields are provided in both the subtool and parent metadata, the metadata in the subtool 
metadata file takes precedence, i.e. subtools inherit metadata from the subtool specific metadata file
and the parent metadata file, with subtool metadata taking precedence.


Tools that are divided into subtools must have a `common/` directory. This directory should contain a parent
metadata file that contains metadata that is pertinent all of the subtools, and may contain requirements files that 
need to be imported by more than one subtool such as SchemaDefRequirement and DockerRequirement.



## Utilities available
Metadata template factory methods

## License
All contributions to this repository shall be made available under the [Apache-2.0](LICENSE.txt). 


## Coming next
- Utilities for generating and validating content. 
- Automated tests.
- Directories and/or repositories for workflows, scripts (another type of command line tool that we'll handle a little 
differently), tool and workflow job files, and input/output file and directory descriptions.


# Todo instances, scripts, workflows.