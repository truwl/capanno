
from ruamel.yaml import safe_load, YAML
from utilities.config import CWL_TOOL_DIR, BASE_DIR, CWL_SCRIPT_DIR, CONTENT_MAPS_DIR
from utilities.classes.script_metadata import ScriptMetadata
from utilities.classes.tool_metadata import ToolMetadata, ParentToolMetadata, SubtoolMetadata
from utilities.helpers.get_paths import get_cwl_tool, get_tool_inputs, get_cwl_tool_metadata, get_tool_version_dir, get_script_version_dir, get_metadata_path, get_relative_path

def make_tools_map(outfile_name='tool-maps'):
    content_map = {}
    for tool_dir in  CWL_TOOL_DIR.iterdir():
        for version_dir in tool_dir.iterdir():
            tool_map = make_tool_map(tool_dir.name, version_dir.name)
            content_map[f"{tool_dir.name}__{version_dir.name}"] = tool_map
    outfile_path = CONTENT_MAPS_DIR / f"{outfile_name}.yaml"
    yaml = YAML(pure=True)
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    with outfile_path.open('w') as outfile:
        yaml.dump(content_map, outfile)
    return


def make_tool_map(tool_name, tool_version):
    tool_map = []
    tool_version_dir = get_tool_version_dir(tool_name, tool_version)
    subdir_names = [subdir.name for subdir in tool_version_dir.iterdir()]
    has_common_dir = True if 'common' in subdir_names else False  # This is the sign of a complex tool. Could also choose len(subdir_names > 1)
    if has_common_dir:
        parent_metadata_path = get_cwl_tool_metadata(tool_name, tool_version, parent=True)
        parent_metadata = ParentToolMetadata.load_from_file(parent_metadata_path)
        parent_rel_path = get_relative_path(parent_metadata_path)
        tool_map.append({'identifier': parent_metadata.identifier, 'path': str(parent_rel_path), 'version': parent_metadata.version, 'name': parent_metadata.name, 'softwareVersion': parent_metadata.softwareVersion, 'type': 'parent'})
        # tool_map[parent_identifier] = str(parent_metadata_path)
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
            subtool_cwl_path = get_cwl_tool(tool_name, tool_version, subtool_name=subtool_name)
            subtool_rel_path = get_relative_path(subtool_cwl_path)
            subtool_metadata_path = get_cwl_tool_metadata(tool_name, tool_version, subtool_name=subtool_name, parent=False)
            subtool_metadata = SubtoolMetadata.load_from_file(subtool_metadata_path)
            tool_map.append({'identifier': subtool_metadata.identifier, 'path': str(subtool_rel_path), 'name': subtool_metadata.name, 'version': subtool_metadata.version, 'type': 'subtool'})
    else: # Not a complex tool. Should just have one directory for main tool.
        metadata_path = get_cwl_tool_metadata(tool_name, tool_version)
        metadata = ToolMetadata.load_from_file(metadata_path)
        tool_rel_path = get_relative_path((get_cwl_tool(tool_name, tool_version)))
        tool_map.append({'identifier': metadata.identifier, 'path': str(tool_rel_path), 'name': metadata.name, 'softwareVersion': metadata.softwareVersion, 'version': metadata.version, 'type': 'tool'})

    return tool_map


def make_script_maps(outfile_name="script-maps"):
    script_maps = {}
    for group_dir in CWL_SCRIPT_DIR.iterdir():
        for project_dir in group_dir.iterdir():
            for version_dir in project_dir.iterdir():
                script_maps[f"{group_dir.name}_{project_dir.name}_{version_dir.name}"] = make_script_map(group_dir.name, project_dir.name, version_dir.name)
    outfile_path = CONTENT_MAPS_DIR / f"{outfile_name}.yaml"
    yaml = YAML(pure=True)
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    with outfile_path.open('w') as outfile:
        yaml.dump(script_maps, outfile)
    return


def make_script_map(group_name, project_name, version):
    script_map = {}
    script_ver_dir = get_script_version_dir(group_name, project_name, version)
    for script_dir in script_ver_dir.iterdir():
        if script_dir.name == 'common':
            continue
        script_cwl_path = script_dir / f"{script_dir.name}.cwl"
        metadata_path = get_metadata_path(script_cwl_path)
        script_metadata = ScriptMetadata.load_from_file(metadata_path)
        script_map[script_metadata.identifier] = {'path': str(script_cwl_path), 'name': script_metadata.name, 'softwareVersion': script_metadata.softwareVersion, 'version': script_metadata.version}
    return script_map



def combine_yaml_files_into_dict(file_path, *file_paths):
    # make a 'master map' that contains all file
    combined_dict = {}
    with file_path.open('r') as f:
        combined_dict.update(safe_load(f))
    for file in file_paths:
        with file.open('r') as f:
            combined_dict.update(safe_load(f))
    return combined_dict

def make_master_map(file_name, *file_names, outfile_name="master_map"):

    file_path = CONTENT_MAPS_DIR / f"{file_name}.yaml"

    raise NotImplementedError

def get_tool_map():
    raise NotImplementedError

def add_to_map(path):
    """
    Add tools scripts and workflows with versions >= 1.0 to content_maps tool_maps, script_maps, or 'workflow_maps'.
    These maps provide an accessible way to work with content based on their identifiers.
    :param path:
    :return:
    """
    raise NotImplementedError