## Repository scope and purpose

This repository contains [Common Workflow Language](https://www.commonwl.org/) (CWL) files, cwl job input files, and associated metadata to describe and run command line tools and workflows with a an emphasis on bioinformatics methods. 

The goal is to create a set CWL tool wrappers that are thorough, well documented, consistent, openly curated, and widely used by the community. A quality set of standardized CWL CommandLineTools will allow for assembly of workflows from "known-good" components with predictable input and output names and behaviour.

CWL files in this repository are divided broadly into tools, scripts, and workflows; collectively referred to as 'methods'. Tools and scripts are both instances of the CommandLineTool CWL class, but differ in the format of their metadata and the software that they describe. Tools typically describe stand-alone programs and are more likely to have documentation and versioned releases. Scripts will often call tools, are more specialized, and have parent and child relationships with other scripts. Workflows correspond to CWL workflows and call one or more tools and/or scripts.

Each CWL/metadata file combination is meant to describe a single version of a method and is separated into its own directory as described in [Repository structure and directory/file names](#structure) below. There are 'instances' directories in each method directory where inputs files (i.e. job files) and metadata that describe individual runs are placed.

### Integration with truwl.com

Methods contributed to this repository that meet a minimum version requirement will be made available through [truwl.com](https://truwl.com). Each method will be viewable from its own web page and be made more findable and accessible through search engine queries. Web pages include comments sections and allow logged in users to save their favorite methods and vote on comments. Each method and use case (job/instance) will also be assigned a unique identifer that can be explicitly referenced.  Relationships between, tools, scripts, workflows, and thier input/output files can be explored and CWL files can be downloaded directly from the site. Workflows will be displayed as interactive graphs from which component tools, scripts, subworkflows, and inputs/outputs can be explored.

## How to contribute

Contributions to CWL files, metadata files, best practices and other documentation files are highly encouraged and can be made through pull requests and issues. CWL and metadata file contributions must follow the best practices described in  [best practices](docs/templates/CommandLineTool_guide.md). 

## <a name="structure"></a> Repository structure and directory/file names

This repository has three main content directories: [cwl-tools](#cwl-tools), [cwl-scripts](#cwl-scripts), and [cwl-workflows](#cwl-workflows).

### <a name="cwl-tools"></a> cwl-tools directory structure

The cwl-tools directory structure follows the pattern `cwl-tools/{toolName}/{toolVersion}/...`

The `{toolName}` path component must be named with the name of the tool being described preserving all capitalization, dashes, and underscores used in the tool's name.

The `{toolVersion}` path component must be named according to the version of the tool that is described preserving the version conventions of the tool author. e.g. v3.2, 3.2, 3b, etc (see [ToDo](docs/components...) for more details).  Patch versions may be omitted. Separate directories should be created for each version of the tool even if the contents of the directories are initially identical to allow the contents to diverge.

####  <a name="tool-subdirectories"></a> Subdirectories of  `cwl-tools/{toolName}/{toolVersion}/`

Tools that contain a limited number of arguments that can be provided together should be described by a single cwl file and metadata file that are contained in a directory named with the tool name. The CWL file must be named `{toolName}.cwl` and the metadata file must be named `{toolName}-metadata.yaml`. e.g. `cwl-tools/cat/8.25/cat/cat.cwl` and `cwl-tools/cat/8.25/cat/cat-metadata.yaml` 

Tools that contain multiple subtools or modes should be divided into logical components (subtools) and be placed in their own directory named `{toolName}-{subtoolName}`e.g. `cwl-tools/STAR/2.7/STAR-genomeGenerate/STAR-genomeGenerate.cwl` and `cwl-tools/STAR/2.7/STAR-genomeGenerate/STAR-genomeGenerate-metadata.yaml`If dividing subtools into more finite components is deemed necessary, more directory nesting must not be added. The addtional subdivision must be specified in the filename e.g. `cwl-tools/{tool_name}/{tool_version}/{tool_name}-{subtool_name}/{tool_name}-{subtool_name}-{subtool_component_name}.cwl`

The metadata file for a subtool should be specific to the subtool. Metadata pertinent to all subtools should be placed in a parent metadata file in the  the `common/` directory described below. If particular metadata fields are provided in both the subtool and parent metadata, the metadata in the subtool metadata file takes precedence, i.e. subtools first get metadata from the parent metadata file then the subtool specific metadata file with the subtool specific metadata overwriting any fields that are provided in both.


Tools that are divided into subtools must have a `common/` directory. This directory should contain a parent metadata file that contains metadata that is pertinent all of the subtools, and may contain requirements files that need to be imported by more than one subtool such as SchemaDefRequirement and DockerRequirement.

##### <a name="tool-instances"></a> Tool instances directories

Each `cwl-tools/{toolName}/{toolVersion}/{tool or subtool name}` directory may have an `instances` directory. This directory will contain any inputs/job files for the tool or subtool and associated metadata about the run. The names of these files should be a four character hexadecimal string  which contributors may provide, or the (arbitrary) filename will be updated by the repository owners to comply with this format. Instance metadata files must have the same name as the job file with -metadata appendend to the base name: e.g. `cwl-tools/cat/8.25/cat/instances/8a6c.yaml` and `cwl-tools/cat/8.25/cat/instances/8a6c-metadata.yaml`. To generate a four character hexadecimal string using Python 3:

~~~Python3
import uuid
uuid.uuid4().hex[:4]
~~~


### <a name="cwl-scripts"></a> cwl-scripts

The cwl-scripts directory structure follows the pattern `cwl-scripts/{groupName}/{projectName}/{version}/{scriptName}/...`

The `groupName` and `projectName` path components should be meaningful and unique enough to not clash with other names.

The `{version}` path component must be named according to the version of the script or group of scripts that are described preserving the version conventions of the script author. (see [ToDo](docs/components...) for more details).

Several `{scriptName}` directories may be present in the `{version}` directory. Generally, scripts that are in the same code repository and that are versioned together should be in the same `{version}` directory. When there are multiple scripts grouped together, there should also be a `cwl-scripts/{groupName}/{projectName}/{version}/common/` directory to place metadata that is relevant to multiple scripts in the directory. See `cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1` for an example.

####  Subdirectories of  `cwl-scripts/{groupName}/{projectName}/{version}/...`

These subdirectories follow the same pattern as [tools subdirectories](#tool-subdirectories)

##### <a name="script-instances"></a> Script instances directories

This follows the same pattern as [tool instances directories](#tool-instances)
 

### <a name="cwl-workflows"></a> cwl-workflows

The `cwl-workflows` directory has a similar structure to `cwl-scripts` and follows the pattern `cwl-workflows/{groupName}/{projectName}/{version}/...`.

### Contents of `cwl-workflows/{groupName}/{projectName}/{version}/...`

Contains workflow CWL file, workflow metadata file, and instances directory.


## Utilities available

Code for generating, validating, and working with files (including generation of tool metadata from [bio.tools](https://bio.tools/)) is in the utilities directory have been developed and tested using Python 3.6 (we like f-strings and use them liberally!). Use of these utilities and setting up a virtual environment is covered in the [Adding Content]() documentation.

## License

All contributions to this repository shall be made available under the [Apache-2.0](LICENSE.txt). 


## What we're working on now (and you could help!)
- More utilities for generating and validating content (CWL, job files, job metadata, etc.)
- Continuous integration
- Referencing software containers.
- utilities for generating CWL tool wrappers.
- More content!

## Code of Conduct

Users of this repo must adhere to xD Bio Inc's [Code of Conduct](https://truwl.com/conduct).