
from pathlib import Path

def get_inputs_schema_template():

    schema_template_path = Path.cwd() / 'tests/test_files/schema-salad/inputs_schema_template.yml'

    return schema_template_path

def get_cwl_tool(tool_name, tool_version, subtool_name=None):
    cwl_tool_path = Path.cwd() / 'cwl-tools' / tool_name / tool_version

    if subtool_name:
        cwl_tool_path = cwl_tool_path / f"{tool_name}_{subtool_name}" / f"{tool_name}-{subtool_name}.cwl"
    else:
        cwl_tool_path = cwl_tool_path / f"{tool_name}" / f"{tool_name}.cwl"

    return cwl_tool_path

def get_cwl_script():
    raise NotImplementedError

def get_tool_inputs(tool_name, tool_version, input_hash, subtool_name=None):

    cwl_tool_dir = get_cwl_tool(tool_name,tool_version, subtool_name=subtool_name).parent
    inputs_path = cwl_tool_dir / 'instances' / f"{input_hash}.yaml"

    return inputs_path


def get_script_inputs():
    raise NotImplementedError