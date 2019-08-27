
from pathlib import Path
from urllib.parse import urlparse
from pprint import pprint
from ruamel.yaml import YAML, safe_load
from schema_salad.main import main as schema_salad_tool
import schema_salad
from cwltool.process import shortname
from utilities.classes.cwl.command_line_tool import load_document
from utilities.helpers.dict_tools import get_dict_from_list
from utilities.helpers.get_files import get_inputs_schema_template, get_cwl_tool, get_tool_inputs, get_cwl_script, get_script_inputs


def validate_tool_inputs(tool_name, tool_version, input_hash, subtool_name=None):

    cwl_file_path = get_cwl_tool(tool_name, tool_version, subtool_name=subtool_name)
    inputs_schema_template = get_inputs_schema_template()
    inputs = get_tool_inputs(tool_name, tool_version, input_hash, subtool_name=subtool_name)

    cwl_document = load_document(str(cwl_file_path))

    cwl_inputs = cwl_document.inputs

    schema_def_requirement = cwl_document._get_schema_def_requirement()

    fields_for_inputs_schema = make_input_fields_for_schema(cwl_inputs, schema_def_requirement)
    with inputs_schema_template.open('r') as template:
        template_dict = safe_load(template)
        _, inputs_field_index = get_dict_from_list(template_dict['$graph'], 'name', 'InputsField')
        template_dict['$graph'][inputs_field_index]['fields'] = fields_for_inputs_schema

    test_file = Path.cwd() / 'temp' / 'test.yaml'

    yaml = YAML(pure=True)
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)

    with test_file.open('w') as tf:
        yaml.dump(template_dict, tf)

    schema_salad_validate(test_file, inputs)
    # return_code = schema_salad_tool(['--debug', str(test_file), str(inputs)])
    # if return_code != 0:
    #     raise ValueError("Placeholder error")
    # print(f"return code: {return_code}")



    return



def make_input_fields_for_schema(cwl_inputs, schema_def_requirement):

    inputs_fields = {}
    for input in cwl_inputs:
        inputs_fields[shortname(input.id)] = {'type': None}
        inputs_fields[shortname(input.id)]['type'] = input._handle_input_type_field(input.type, schema_def_requirement)
    return inputs_fields


def schema_salad_validate(schema_path, document_path):
    '''
    Adapted from schema_salad main().
    :param schema_path:
    :param document_path:
    :return:
    '''

    strict_foreign_properties = False
    strict = True
    metaschema_names, metaschema_doc, metaschema_loader = schema_salad.main.schema.get_metaschema()
    schema_uri = str(schema_path)
    if not (urlparse(schema_uri)[0] and urlparse(schema_uri)[0] in ['http', 'https', 'file']):
        schema_uri = schema_salad.main.file_uri(schema_uri)
    schema_raw_doc = metaschema_loader.fetch(schema_uri)

    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, schema_uri)

    # Validate schema against metaschema
    schema_salad.main.schema.validate_doc(metaschema_names, schema_doc, metaschema_loader, True)

    # Get the json-ld context and RDFS representation from the schema
    metactx = schema_salad.main.schema.collect_namespaces(schema_metadata)
    if "$base" in schema_metadata:
        metactx["@base"] = schema_metadata["$base"]

    (schema_ctx, rdfs) = schema_salad.main.jsonld_context.salad_to_jsonld_context(
        schema_doc, metactx)

    # Create the loader that will be used to load the target document.
    document_loader = schema_salad.main.Loader(schema_ctx, skip_schemas=False)

    # Make the Avro validation that will be used to validate the target
    # document

    avsc_obj = schema_salad.main.schema.make_avro(schema_doc, document_loader)

    avsc_names = schema_salad.main.schema.make_avro_schema_from_avro(avsc_obj)

    # Load target document and resolve refs
    uri = str(document_path)
    document, doc_metadata = document_loader.resolve_ref(uri, strict_foreign_properties=strict_foreign_properties, checklinks=False) # This is what's getting us around file link checking.

    schema_salad.main.schema.validate_doc(avsc_names, document, document_loader, strict=strict, strict_foreign_properties=strict_foreign_properties)

    print(f"{document_path} is valid.")
    return