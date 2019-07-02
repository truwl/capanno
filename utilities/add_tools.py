#!/usr/bin/env python3

from pathlib import Path
from utilities.classes.tool_metadata import ToolMetadata, ParentToolMetadata


def add_parent_tool(tool_name, tool_version, subtool_names=None, biotools_id=None):
    """
    Make the correct directory structure for adding a new command line tool. Optionally, create initialized CWL
    and metadata files. Run from cwl-tools directory.
    :param tool_name(str): Name of the tool
    :param tool_version(str): version of the tool
    :param subtool_names(list(str)): list of subtool names if the tool is broken into multiple subtools.
    :param mk_meta_files(bool): Specify whether to make initial CWL and metadata files.
    :return: None
    """
    tool_version = str(tool_version) # In case ArgumentParser is bypassed.

    common_dir = Path().cwd() / 'cwl-tools'/ tool_name / tool_version / 'common'
    common_dir.mkdir(parents=True)
    if biotools_id:
        parent_metadata = ParentToolMetadata.create_from_biotools(biotools_id, tool_version, subtool_names)
    else:
        parent_metadata = ParentToolMetadata(name=tool_name, softwareVersion=tool_version, featureList=subtool_names)
    new_file_path = common_dir / f"{parent_metadata.name}-metadata.yaml"
    parent_metadata.mk_file(new_file_path)
    return new_file_path


def add_tool(tool_name, tool_version, biotools_id=None):
    tool_version = str(tool_version)  # In case ArgumentParser is bypassed.
    new_tool_dir = Path.cwd() / 'cwl-tools' / tool_name / tool_version / tool_name
    new_tool_dir.mkdir(parents=True)
    if biotools_id:
        new_tool = ToolMetadata.create_from_biotools(biotools_id, tool_version)
    else:
        new_tool = ToolMetadata(name=tool_name, softwareVersion=tool_version)

    new_tool_path = new_tool_dir / f"{tool_name}-metadata.yaml"
    new_tool.mk_file(new_tool_path)
    return new_tool_path


def add_subtool(parent_rel_path, subtool_name):
    parent_path = Path(parent_rel_path) # Should be in /common
    parent_meta = ParentToolMetadata.load_from_file(parent_path)
    new_subtool_dir = parent_path.parents[1] / f"{parent_meta.name}_{subtool_name}"
    new_subtool_dir.mkdir()
    subtool_meta = parent_meta.make_subtool_metadata(subtool_name)
    new_file_path = new_subtool_dir / f"{subtool_name}-metadata.yaml"
    subtool_meta.mk_file(new_file_path)
    return new_file_path
