# Getting Started

## General

THIS IS OUT OF DATE. Tools to manage content have been factored out and released on pypi as xd-cwl-utils !!!!! Documentation on using that tool will come soon (soon is a relative term.).

To start working with this repo you should first clone it and create a python 3.6+ virtual environment.

~~~bash
$ git clone https://github.com/xDBio-Inc/cwl-tools.git
$ cd cwl-tools
~~~

Create a virtual environment with virtualenvwrapper and install requirements.
~~~bash
$  mkvirtualenv -p <path to python 3.6+> <virtual env name>  # specify python version and name of new environment
(venv name) $ pip install -r requirements.txt
~~~

## CWL
CWL files must follow [best practices]() described in this repository and be validated using [cwltool](https://github.com/common-workflow-language/cwltool). CWL validation may soon be incorporated into this repository.

## Metadata

### Adding tool metadata
To add a new standalone tool (metadata file and directory structure) from terminal:
~~~bash
$ python -m utilities.add_content tool <tool name> <tool version>
~~~

To add a new standalone tool from terminal and initialize metadata from [bio.tools](https://bio.tools):

~~~bash
$ python -m utilities.add_content tool <tool name> <tool version> --biotoolsID <biotools id>
~~~

To add new parent tool metadata

~~~bash
$ python -m utilities.add_content parent_tool <tool name> <tool version> <subtool names> [--biotoolsID <biotools id>]
~~~

To add new subtool

~~~bash
$ python -m utilities.add_content subtool <path to parent metadata> <subtool name>
~~~

Tools can also be added from the python REPL.

~~~python3
>>> from utilities.add.add_tools import add_tool
>>> add_tool(tool_name, tool_version, biotools_id='star')
~~~


### Updating tool metadata

The simplest way to update tool metadata is to edit the metadata files, then validate that it is formatted correctly.


### Validate tool metadata

To validate a metadata file from the command line

~~~bash
$ python -m utilities.validate [tool | parent_tool | subtool] <path to metadata file>
~~~


## Scripts
To add a new script (metadata file and directory structure)

~~~bash
$ python -m utilities.add_content script <group name> <project name> <script_version> <script_name>
~~~

To add new common script metadata (metadata that is common to multiple scripts)
~~~bash
$ python -m utilities.add_content common_script <group name> <project name> <script_version> <filename>
~~~

### Updating script metadata

The simplest way to update script metadata is to edit the metadata files, then validate that it is formatted correctly.


### Validate script

To validate a metadata file from the command line

~~~bash
$ python -m utilities.validate [script | common_script] <path to metadata file>
~~~

## Workflows

### Adding workflow metadata

To add a new workflow (metadata file and directory structure)

~~~bash
$ python -m utilities.add_content workflow <group name> <workflow name> <workflow version>
~~~

### Updating workflow metadata

The simplest way to update workflow metadata is to edit the metadata files, then validate that it is formatted correctly.

### Validate workflow

To validate a metadata file from the command line

~~~bash
$ python -m utilities.validate workflow <path to metadata file>
~~~