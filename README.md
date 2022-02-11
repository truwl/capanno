[![Build Status](https://travis-ci.com/truwl/capanno.svg?branch=master)](https://travis-ci.com/truwl/capanno)
## Repository scope and purpose

This repository contains metadata, workflow language files, and corresponding inputs files that describe jobs for computational tools, scripts, and workflows--collectively referred to as 'methods'--used in bioinformatics. The repository focuses on four major pipeline frameworks: [Common Workflow Language](https://www.commonwl.org/) (CWL), [Workflow Description Language](https://openwdl.org/) (WDL), [Nextflow](http://nextflow.io/), and [Snakemake](https://snakemake.github.io/).

The goal is to document, collect resources, and openly curate methods and usage examples to facilitate easier access and use by life science researchers. 

Resources in this repository are divided broadly into three directories: tools, scripts, and workflows. Tools and scripts differ in the format of their metadata and the software that they describe. Tools typically consist stand-alone programs and are more likely to have documentation and versioned releases. Scripts will often call tools, are more specialized, and have parent and child relationships with other scripts--they call or are called by other scripts. Workflows correspond to series of tools and/or scripts that are called and executed in an interdependent manner.

Each method directory has an `instances` directory where inputs files that describe a job and metadata that describe individual runs can be placed.

### Integration with truwl.com
Validated methods contributed to this repository will be made available through [truwl.com](https://truwl.com). Each method will be viewable from its own web page and be made more findable and accessible through search engine queries. Web pages include comments sections and allow logged in users to save their favorite methods and vote on comments. Each method and usage example (job/instance) will also be assigned a unique identifier that can be explicitly referenced. Relationships between, tools, scripts, workflows, and thier input/output files can be explored and tool wrapper files can be downloaded directly from the site. Some workflows can be displayed as interactive graphs from which component tools, scripts, subworkflows, and inputs/outputs can be explored. Many of the workflows described in WDL can be executed, monitored, and shared directly from the site. 

## Repository management with `capanno-utils`
[`capanno-utils`](https://github.com/truwl/capanno-utils) is a companion Python package to simplify generating and validating content in this repository. Basic use of `capanno-utils` is described in [Getting Started](docs/Getting_Started.md)

## How to contribute

Contributions including additions to methods as well as corrections and enhancements to existing files and documentation are highly encouraged and can be made through pull requests and issues. See [Getting Started](docs/Getting_Started.md) to learn how to add new content with the proper directory structure and initialize metadata files. 

## <a name="structure"></a> Repository structure and directory/file names

This repository has three main content directories: [tools](#tools), [scripts](#scripts), and [workflows](#workflows).
```
capanno/
├── scripts
├── tools
└── workflows
```

### <a name="tools"></a> tools directory structure

The directory structure for a tool is
```
tools/{toolName}
├── {versionName 1}
│   ├── common
│   │   ├── common-metadata.yaml
│   │   ├── common_file_to_import.yaml
│   │   └── other_common_file_to_import.yml
│   ├── toolName
│   │   ├── instances/
│   │   ├── toolName-subtoolName1-metadata.yaml
│   │   └── toolName-subtoolName1.[cwl | wdl | nf | snakefile]
│   ├── toolName-subtoolName1
│   │   ├── instances/
│   │   ├── toolName-subtoolName1-metadata.yaml
│   │   └── toolName-subtoolName1.[cwl | wdl | nf | snakefile]
│   └── toolName_subtoolName2
│       ├── instances/
│       ├── toolName-subtoolName2-metadata.yaml
│       └── toolName-subtoolName2.[cwl | wdl | nf | snakefile]
└── versionName 2/...

```

This structure and initialized metadata files are easily generated using `capanno-utils` when adding a new tool, versionName directory, or subtool.

The `{toolName}` path component must be named with the name of the tool being described preserving all capitalization, dashes, and underscores used in the tool's name.

The `{versionName}` path component should be named according to version conventions of the tool author as much as possible. Since metadata and workflow language files can be applicable to multiple versions of the tool, version names can, and usually do contain variable portions; e.g. 1.x, 1.2.x, etc. A list of specific versions that are covered by the version name can be specified in the `common-metadata.yaml` file. When tool descriptions, other metadata, and parameters can no longer accurately be described by a single set of files a new `{versionName}` directory should be created to allow the differences in versions to diverge.

The scope of functions that a tool performs can vary widely. Tools may be developed as a single program or may include multiple subtools or modes that are specified when called, e.g. `toolName subtoolName ...`. Metadata and relevant inputs and parameters differ depending on the subtool. Tools that contain multiple subtools or modes are divided into these logical components and described in their own subdirectories named `{toolName}-{subtoolName}` to accommodate these differences e.g. `tools/STAR/2.7/STAR-genomeGenerate/`. For tools that can be called without specifying a subtool, a {toolName} directory is included e.g `tools/cat/8.x/cat/`. If dividing subtools into more finite components is necessary, additional subdivision must be specified in the filenames e.g. `tools/{tool_name}/{tool_version}/{tool_name}-{subtool_name}/{tool_name}-{subtool_name}-{subtool_component_name}-metadata.yaml` rather than adding more directory nesting. Each tool/version combination contains a `common/` directory. This directory must contain a `common-metadata.yaml` file that contains metadata that is pertinent all of the subtools, and may contain additional files that are pertinent to multiple subtools such as requirements files that need to be imported by more than one subtool such as SchemaDefRequirement and DockerRequirement descriptions for CWL files.

If the same metadata fields are provided in both the subtool and common metadata, the metadata in the subtool metadata file takes precedence.

##### tool metadata

The metadata fields for tools are described in the [tool metadata guide](docs/ToolGuide.md)

##### <a name="tool-instances"></a> Tool instances directories

Each `tools/{toolName}/{toolVersion}/{tool or subtool name}` directory has an `instances` directory. This directory will contain any inputs files for the tool or subtool and associated metadata about the run. Instance metadata files must have the same name as the job file with -metadata appended to the base name: e.g. `tools/cat/8.25/cat/instances/8a6c.yaml` (job file) and `tools/cat/8.25/cat/instances/8a6c-metadata.yaml` (metadata for job file). 


### <a name="scripts"></a> scripts

The directory structure for scripts is arranged around group and project names, similar to GitHub repositories.

```
scripts/{groupName}/
├── {projectName1}
│   └── {versionName1}
│       └── {scriptName1}
│           ├── {scriptName1}-metadata.yaml
│           ├── {scriptName1}.[cwl | wdl | nf | snakefile]
│           └── instances/
└── {projectName2}
    ├── {versionName2}
    │   ├── common
    │   │   └── common-metadata.yaml
    │   ├── {scriptName2}
    │   │   ├── {scriptName2}-metadata.yaml
    │   │   ├── {scriptName2}.[cwl | wdl | nf | snakefile]
    │   │   └── instances/
    │   └── {scriptName3}
    │   │   ├── {scriptName3}-metadata.yaml
    │   │   ├── {scriptName3}.[cwl | wdl | nf | snakefile]
    │       └── instances/
    └── {versionName3}
        ├── common
        │   └── common-metadata.yaml
        ├── {scriptName2}
        │   ├── {scriptName2}-metadata.yaml
        │   ├── {scriptName2}.[cwl | wdl | nf | snakefile]
        │   └── instances/
        └── {scriptName3}
            ├── {scriptName3}-metadata.yaml
            ├── {scriptName3}.[cwl | wdl | nf | snakefile]
            └── instances/
```

The `groupName` and `projectName` path components should be meaningful and unique enough to not clash with other names. Each group can have many projects.

The `{versionName}` path component must be named according to the version of the script or group of scripts that are described preserving the version conventions of the script author.

Several `{scriptName}` directories may be present in each `{versionName}` directory. Generally, scripts that are in the same code repository and that are versioned together should be in the same `{versionName}` directory. When there are multiple scripts grouped together, there should also be a `scripts/{groupName}/{projectName}/{version}/common/` directory to place metadata that is relevant to multiple scripts in the directory. 

##### Script metadata
The metadata fields for tools are described in the [script metadata guide](docs/ScriptGuide.md)

##### <a name="script-instances"></a> Script instances directories

This follows the same pattern as [tool instances directories](#tool-instances)


### <a name="workflows"></a> workflows

The `workflows` directory has a similar structure to `scripts` and is arranged around group and project names.

Contains workflow file, workflow metadata file, and instances directory.

```
workflows/{groupName}/
├── {projectName1}
│   ├── {versionName1a}
│   │   ├── {empty inputs specifation file}
│   │   ├── {projectName1}-metadata.yaml
│   │   ├── {workflowName}.[cwl | wdl | nf | snakefile]
│   │   └── instances
│   │       ├── {instanceName}-metadata.yaml
│   │       └── {instanceName}.json
│   ├── versionName1b
│   │   ├── {empty inputs specifation file}
│   │   ├── {projectName}-metadata.yaml
│   │   ├── {workflowName}.[cwl | wdl | nf | snakefile]
│   │   └── instances/
├── {projectName2}
│   └── versionName2a
│       ├── {projectName2}-metadata.yaml
│       ├── additional-files...
│       └── instances/
```

##### Workflow metadata
The metadata fields for workflows are described in the [workflow metadata guide](docs/WorkflowGuide.md)

## License

All contributions to this repository shall be made available under the [Apache-2.0](LICENSE.txt). 


## What we're working on now (and you could help!)
- More utilities for generating and validating content (job files, job metadata, etc.)
- Referencing software containers.
- utilities for generating tool wrappers.
- More content!

## Code of Conduct

Users of this repo must adhere to Truwl's [Code of Conduct](https://github.com/xDBio-Inc/Policies/blob/master/xd_code_of_conduct.md).
