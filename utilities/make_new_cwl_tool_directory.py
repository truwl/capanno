#!/usr/bin/env python3

from os import makedirs, mkdir
import sys
import argparse
from pathlib import Path
from os.path import join as pjoin
from shutil import copyfile
from utilities.classes.metadata import ToolMetadata, SubtoolMetadata, ParentToolMetadata
from utilities.get_metadata import make_metadata_file_from_biotools


parser = argparse.ArgumentParser(description='Initialize directory and files for new cwl-tool')
parser.add_argument('name', type=str, nargs=1, help='The name of the new tool that is being added')
parser.add_argument('version', type=str, nargs=1, help='The version of the new tool that is being added')
parser.add_argument('subtools', type=str, nargs='*', help='List of subtool names for the tool.')
parser.add_argument('--biotools_id', help="bio.tools ID. If provided will be used to populate metadata file with data from https://bio.tools.")

def mk_tool_directory(tool_name, tool_version, subtool_names=None, biotools_id=None):
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
    # base_path = Path().cwd()
    if subtool_names:
        common_dir = Path.cwd() / tool_name / tool_version / 'common'
        common_dir.mkdir()
        if biotools_id:
            make_metadata_file_from_biotools(biotools_id, common_dir)
        else:
            parent_metadata = ParentToolMetadata(name=tool_name, softwareVersion=tool_version)
        for subtool_name in subtool_names:
            subtool_path = Path().cwd() / 'cwl-tools' / tool_name / tool_version / f"{tool_name}_{subtool_name}"
            subtool_path.mkdir(parents=True)


        if biotools_id:
            subtool_template = template_dir + '/subtool-metadata-template.yaml'
            parent_template = template_dir + '/parent-tool-metadata-template.yaml'
            [copyfile(subtool_template, pjoin(tool_name, tool_version, tool_name + '_' + s_name, tool_name + '-' + s_name + '-metadata.yaml')) for s_name in subtool_names]
            copyfile(parent_template, pjoin(tool_name, tool_version, 'common', tool_name + '-metadata.yaml'))
        if len(subtool_names) == 1:
            mkdir(pjoin(tool_name, tool_version, tool_name))

    else:
        makedirs(pjoin(tool_name, tool_version, tool_name))
        if mk_meta_files:
            tool_template = template_dir + '/tool-metadata-template.yaml'
            copyfile(tool_template, pjoin(tool_name, tool_version, tool_name, tool_name + '-metadata.yaml'))
    return


def _pop_fields(metafile_path, tool_name, subtool_name, tool_version):
    # TODO populate name, softwareVersion, applicationSuite, etc. when iniitalizing metadata files.
    # populate metadata files from mk_tool_directory arguments.
    raise NotImplementedError

if __name__ == "__main__":
    args = parser.parse_args()
    print(args.name, args.version, args.subtools)
    # mk_tool_directory(args.name, args.version, args.subtool_names, args.biotools_id)
    sys.exit(0)