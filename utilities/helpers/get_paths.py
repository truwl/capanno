
from pathlib import Path
from utilities.config import CWL_TOOL_DIR

# Misc

def get_inputs_schema_template():

    schema_template_path = Path.cwd() / 'tests/test_files/schema-salad/inputs_schema_template.yml'

    return schema_template_path

# cwl-tools

def get_tool_version_dir(tool_name, tool_version):
    version_dir = CWL_TOOL_DIR / tool_name / tool_version
    return version_dir


def get_cwl_tool(tool_name, tool_version, subtool_name=None):
    version_dir = get_tool_version_dir(tool_name, tool_version)

    if subtool_name:
        cwl_tool_path = version_dir / f"{tool_name}_{subtool_name}" / f"{tool_name}-{subtool_name}.cwl"
    else:
        cwl_tool_path = version_dir / f"{tool_name}" / f"{tool_name}.cwl"
    return cwl_tool_path


def get_cwl_tool_metadata(tool_name, tool_version, subtool_name=None, parent=False):
    version_dir = get_tool_version_dir(tool_name, tool_version)

    if parent:
        cwl_tool_metadata_path = version_dir / 'common' / f"{tool_name}-metadata.yaml"
    else:
        cwl_tool_metadata_path = get_metadata_path(get_cwl_tool(tool_name, tool_version, subtool_name=subtool_name))
    return cwl_tool_metadata_path

# cwl-scripts

def get_cwl_script():
    raise NotImplementedError

def get_tool_inputs(tool_name, tool_version, input_hash, subtool_name=None):

    cwl_tool_dir = get_cwl_tool(tool_name,tool_version, subtool_name=subtool_name).parent
    inputs_path = cwl_tool_dir / 'instances' / f"{input_hash}.yaml"

    return inputs_path


def get_script_inputs():
    raise NotImplementedError


# cwl-workflows


# helpers

def get_metadata_path(cwl_path):
    path_dir = cwl_path.parent
    metafile_name = f"{cwl_path.stem}-metadata.yaml"
    metadata_path = path_dir / metafile_name
    return metadata_path
