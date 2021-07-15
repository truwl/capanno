# Getting Started

## General

`capanno-utils` is a package for managing content in this repository and is specified in `requirements.txt`

To start working with this repo you should first clone it and create a python 3.6+ virtual environment.

~~~bash
$ git clone https://github.com/truwl/capanno.git
$ cd capanno
~~~

Create a virtual environment, shown here using [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/), and install requirements.
~~~bash
$  mkvirtualenv -p <path to python 3.6+> <venv name>  # specify python version and name of new environment
(venv name) $ pip install -r requirements.txt # or can use $pip install capanno-utils
~~~

##Adding new content
### tools
To add a new tool use the `capanno-utils` commands `capanno-add tool` and `capanno-add subtool`.

~~~shell
usage: capanno-add tool [-h] [--biotoolsID BIOTOOLSID] [--has-primary]
                        [--init-cwl] [--init-wdl] [--init-nf] [--init-sm]
                        [--no-clobber]
                        tool_name version_name
                        [subtool_names [subtool_names ...]]

positional arguments:
  tool_name             The name of the tool.
  version_name          The version name of the tool.
  subtool_names         The names of the subtools/subcommands for the tool.
                        Directories and files will be initialized for the
                        specified subtools.

optional arguments:
  -h, --help            show this help message and exit
  --biotoolsID BIOTOOLSID
                        biotools id from https://bio.tools. Metadata fields
                        will be populated with data from bio.tools if
                        specified
  --has-primary         Specify if tool is callable without a subcommand.
  --init-cwl            If specified, CWL CommandLineTool files will be
                        intiated for the subtools and primary tool if it
                        exists.
  --init-wdl            If specified, wdl task files will be intiated for the
                        subtools and primary tool if it exists.
  --init-nf             If specified, nextflow files will be intiated for the
                        subtools and primary tool if it exists.
  --init-sm             If specified, snakemake files will be intiated for the
                        subtools and primary tool if it exists.
  --no-clobber          Do not overwrite the existing tool_name, version_name
                        combination in the repository, if it already exists.

~~~

To add a subtool entry to an already existing tool entry use `capanno-add subtool`
```shell
usage: capanno-add subtool [-h] [-u] [--init-cwl [INIT_CWL]]
                           [--init-wdl [INIT_WDL]] [--init-nf [INIT_NF]]
                           [--init-sm [INIT_SM]] [--no-clobber]
                           tool_name version_name subtool_name

positional arguments:
  tool_name             The name of the tool.
  version_name          The version name of the tool.
  subtool_name          The name of the subtool to add. The subtool name must
                        be present in the common-metadata featureList field.
                        The featureList field can be updated automatically
                        with the --update-featureList option.

optional arguments:
  -h, --help            show this help message and exit
  -u, --update-featureList
                        Update the featureList of the Application Suite
                        metadata to contain new subtool.
  --init-cwl [INIT_CWL]
                        If specified, CWL CommandLineTool files will be
                        intiated. If a url is provided, the file will be
                        intialized from the url.
  --init-wdl [INIT_WDL]
                        If specified, WDL task file will be intiated . If a
                        url is provided, the file will be intialized from the
                        url.
  --init-nf [INIT_NF]   If specified, a nextflow file will be intiated . If a
                        url is provided, the file will be intialized from the
                        url.
  --init-sm [INIT_SM]   If specified, a snakemake file will be intiated. If a
                        url is provided, the file will be intialized from the
                        url.
  --no-clobber          Do not overwrite the existing tool_name, version_name,
                        subtool name combination in the repository, if it
                        already exists.

```

####Examples
Add `STAR` to the repository with a version name of 2.x. STAR has subtools alignReads, genomeGenerate, inputAlignmentsFromBAM, and liftOver. There is an entry for STAR in bio.tools with an id of 'star'. The genomeGenerate subtool is intentionally excluded in this example and added in the next example.
```bash
$ capanno-add tool STAR 2.x alignReads inputAlignmentsFromBAM liftOVer --biotoolsID star
```
Creates the directory structure below. `common/common-metadata.yaml` is populated with metadata from bio.tools. Metadata can also be added or updated directly in the metadata file if the bio.tools entry is incomplete or there is not a bio.tools entry for the tool. The metadata fields are described in [Tool Guide](ToolGuide.md). 
```shell
tools/STAR/2.x/
├── STAR_alignReads
│   ├── STAR-alignReads-metadata.yaml
│   └── instances/
├── STAR_inputAlignmentsFromBAM
│   ├── STAR-inputAlignmentsFromBAM-metadata.yaml
│   └── instances/
├── STAR_liftOVer
│   ├── STAR-liftOVer-metadata.yaml
│   └── instances/
└── common
    └── common-metadata.yaml
````

Add the subtool genomeGenerate to the STAR 2.x folder and initialize a cwl file from a url.
```bash
$ capanno-add subtool STAR 2.x genomeGenerate --init-cwl https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/STAR/STAR-Index.cwl -u
```
There is now a `genomeGenerate` directory that contains a CWL file. The specified CWL file must be valid CWL. To add draft versions of CWL files that are not valid CWL, we recommend you use `curl` or similar.
```shell
tools/STAR/2.x/
├── STAR_alignReads
│   ├── STAR-alignReads-metadata.yaml
│   └── instances/
├── STAR_genomeGenerate
│   ├── STAR-genomeGenerate-metadata.yaml
│   ├── STAR-genomeGenerate.cwl
│   └── instances/
├── STAR_inputAlignmentsFromBAM
│   ├── STAR-inputAlignmentsFromBAM-metadata.yaml
│   └── instances/
├── STAR_liftOVer
│   ├── STAR-liftOVer-metadata.yaml
│   └── instances/
└── common
    └── common-metadata.yaml
```
### workflows
To add a new workflow use the `capanno-utils` command `capanno-add workflow`

```shell
usage: capanno-add workflow [-h] group_name workflow_name workflow_version

positional arguments:
  group_name        The name of the group directory that the workflow will go
                    into.
  workflow_name     The name of the workflow. Replace spaces with underscores.
  workflow_version  The version of the workflow.

optional arguments:
  -h, --help        show this help message and exit

```
####Examples
Add a workflow for a group `example_group` with a project name `sequencing_workflow` version `1.fake`.

```shell
capanno-add workflow example_group sequencing_workflow 1.fake
```
This command creates the directory structure below
```
workflows/example_group
└── sequencing_workflow
    └── 1.fake
        ├── instances/
        └── sequencing_workflow-metadata.yaml
```
The workflow file(s) and additional related files can then be added directly to the directory.

The metadata fields for workflows are described in [Workflow Guide](WorkflowGuide.md)

## Workflow Language Files
### CWL
CWL files should follow [best practices]() described in this repository and must be validated using capanno-utils as described below.

### Snakemake
No guidance or validation is provided for Snakemake files at this time.

### Nextflow
No guidance or validation is provided for Nextflow files at this time.

### WDL
WDL files must pass a cursory validation step that checks if a WDL file can be successfully loaded.


## Validating content 

To validate content use the `capanno-utils` command `capanno-validate`

~~~shell
usage: capanno-validate [-h] [-p ROOT_PATH] [-q] path

Validate metadata and workflow language files.

positional arguments:
  path                  Provide the path to validate. If a directory is
                        specified, all content in the directory will be
                        validated. If a file is specified, only that file will
                        be validated.

optional arguments:
  -h, --help            show this help message and exit
  -p ROOT_PATH, --root-repo-path ROOT_PATH
                        Specify the root path of your content repo if it is
                        not the current working directory.
  -q, --quiet           Silence messages to stdout

~~~

#### Examples
Validate all the content in the capanno repository.
```shell
capanno-validate . # Assumes capanno is your current working directory.
```
Validate all the content in the tools directory. Can validate scripts/ or workflows/ directories the same way.
```shell
capanno-validate tools/ # Assumes capanno is your current working directory.
```

Validate a single metadata, or workflow (CWL/WDL) file:
```shell
capanno-validate tools/<tool_name>/<version_name>/<subtool_name>/<file to validate>
```
