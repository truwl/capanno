
from ruamel.yaml import safe_load, YAML
from utilities.config import CWL_TOOL_DIR, BASE_DIR
from utilities.helpers.get_paths import get_cwl_tool, get_tool_inputs, get_cwl_tool_metadata, get_tool_version_dir

def make_content_map():
    content_map = {}
    for tool_dir in  CWL_TOOL_DIR.iterdir():
        for version_dir in tool_dir.iterdir():
            tool_map = make_tool_map(tool_dir.name, version_dir.name)
            content_map[f"{tool_dir.name}_{version_dir.name}"] = tool_map
    outfile_path = BASE_DIR / 'temp' / 'content_map_test.yaml'
    yaml = YAML(pure=True)
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    with outfile_path.open('w') as outfile:
        yaml.dump(content_map, outfile)
    return


def make_tool_map(tool_name, tool_version):
    tool_map = {}
    tool_version_dir = get_tool_version_dir(tool_name, tool_version)
    subdir_names = [subdir.name for subdir in tool_version_dir.iterdir()]
    complex_tool = True if 'common' in subdir_names else False  # This is the sign of a complex tool. Could also choose len(subdir_names > 1)
    if complex_tool:
        tool_metadata = get_cwl_tool_metadata(tool_name, tool_version, parent=True)
        with tool_metadata.open('r') as metadata_file:
            parent_metadata = safe_load(metadata_file)
        parent_version = parent_metadata['version']
        parent_identifier = parent_metadata['identifier']
        tool_map[parent_identifier] = str(tool_metadata)
        subdir_names.remove('common')
        for subdir_name in subdir_names:
            subtool_name_parts = subdir_name.split('_')
            subtool_name_parts_len = len(subtool_name_parts)
            if subtool_name_parts_len == 2:  # The common case.
                tool_name_from_dir, subtool_name = subdir_name.split('_')
            elif subtool_name_parts_len == 1: # have a subtool that is the 'main' part of the tool; i.e. not a submodule. e.g. md5sum and md5sum_check
                tool_name_from_dir = subtool_name_parts[0]
                subtool_name = None
            else:
                raise NotImplementedError(
                    f"There are zero or more than one underscore in directory {subdir_name}. Can't handle this yet.")
            assert (tool_name_from_dir == tool_name), f"{tool_name} should be equal to {tool_name_from_dir}."
            subtool_cwl = get_cwl_tool(tool_name, tool_version, subtool_name=subtool_name)
            subtool_metadata = get_cwl_tool_metadata(tool_name, tool_version, subtool_name=subtool_name, parent=False)
            with subtool_metadata.open('r') as subtool_metadata_file:
                subtool_metadata_dict = safe_load(subtool_metadata_file)
            subtool_version = subtool_metadata_dict['version']
            subtool_identifier = subtool_metadata_dict['identifier']

            tool_map[subtool_identifier] = str(subtool_cwl)
    else: # Not a complex tool. Should just have one directory for main tool.
        tool_metadata = get_cwl_tool_metadata(tool_name, tool_version)
        with tool_metadata.open('r') as metadata_file:
            tool_metadata_dict = safe_load(metadata_file)
        tool_identifier = tool_metadata_dict['identifier']
        tool_cwl_version = tool_metadata_dict['version']

        tool_map[tool_identifier] = str(get_cwl_tool(tool_name, tool_version))
    return tool_map





def add_to_map(path):
    """
    Add tools scripts and workflows with versions >= 1.0 to content_maps tool_maps, script_maps, or 'workflow_maps'.
    These maps provide an accessible way to work with content based on their identifiers.
    :param path:
    :return:
    """
    raise NotImplementedError