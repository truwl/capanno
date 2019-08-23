
from pathlib import Path
from pprint import pprint
from cwltool.process import shortname
from utilities.classes.cwl.command_line_tool import load_document
from utilities.helpers.get_files import get_inputs_schema_template, get_cwl_tool, get_tool_inputs, get_cwl_script, get_script_inputs


def validate_tool_inputs(tool_name, tool_version, input_hash, subtool_name=None):

    cwl_file_path = get_cwl_tool(tool_name, tool_version, subtool_name=subtool_name)
    inputs_schema_template = get_inputs_schema_template()
    inputs = get_tool_inputs(tool_name, tool_version, input_hash, subtool_name=subtool_name)

    cwl_document = load_document(str(cwl_file_path))

    cwl_inputs = cwl_document.inputs

    schema_def_requirement = cwl_document._get_schema_def_requirement()

    fields_for_inputs_schema = make_input_fields_for_schema(cwl_inputs, schema_def_requirement)

    pprint(fields_for_inputs_schema)

    return




def make_input_fields_for_schema(cwl_inputs, schema_def_requirement):

    inputs_fields = {}
    for input in cwl_inputs:
        inputs_fields[shortname(input.id)] = input._handle_input_type_field(input.type, schema_def_requirement)
    return inputs_fields

