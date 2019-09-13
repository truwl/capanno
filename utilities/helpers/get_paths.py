
from pathlib import Path
from utilities.config import CWL_TOOL_DIR, CWL_SCRIPT_DIR, BASE_DIR

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

def get_tool_inputs(tool_name, tool_version, input_hash, subtool_name=None):

    cwl_tool_dir = get_cwl_tool(tool_name,tool_version, subtool_name=subtool_name).parent
    inputs_path = cwl_tool_dir / 'instances' / f"{input_hash}.yaml"

    return inputs_path

# cwl-scripts

def get_script_version_dir(group_name, project_name, version):
    script_ver_dir = CWL_SCRIPT_DIR / group_name / project_name / version
    return script_ver_dir

def get_cwl_script(group_name, project_name, version, script_name):
    script_ver_dir = get_script_version_dir(group_name, project_name, version)
    script_path = script_ver_dir / script_name / f"{script_name}.cwl"
    return script_path

def get_script_metadata(group_name, project_name, version, script_name):
    script_ver_dir = get_script_version_dir(group_name, project_name, version)
    script_metadata_path = script_ver_dir / script_name / f"{script_name}-metadata.cwl"
    return script_metadata_path


def get_script_inputs():
    raise NotImplementedError


# cwl-workflows


# helpers

def get_relative_path(full_path, base_path=BASE_DIR):

    return full_path.relative_to(base_path)

def get_metadata_path(cwl_path):
    path_dir = cwl_path.parent
    metafile_name = f"{cwl_path.stem}-metadata.yaml"
    metadata_path = path_dir / metafile_name
    return metadata_path
