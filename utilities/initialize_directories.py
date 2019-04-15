

from os import makedirs, mkdir, environ
from os.path import join as pjoin
from shutil import copyfile

def mk_tool_directory(tool_name, tool_version, subtool_names=None, mk_meta_files=True):
    """
    Make the correct directory structure for adding a new command line tool. Optionally, create initialized CWL
    and metadata files.
    :param tool_name(str): Name of the tool
    :param tool_version(str): version of the tool
    :param subtool_names(list(str)): list of subtool names if the tool is broken into multiple subtools.
    :param mk_meta_files(bool): Specify whether to make initial CWL and metadata files.
    :return: None
    """
    tool_version = str(tool_version)
    template_dir = environ['PWD'] + '/templates'
    if subtool_names:
        [makedirs(pjoin(tool_name, tool_version, tool_name + '_' + s_name)) for s_name in subtool_names]
        mkdir(pjoin(tool_name, tool_version, 'common'))
        if mk_meta_files:
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


def pop_fields(metafile_path, tool_name, subtool_name, tool_version):
    # TODO populate name, softwareVersion, applicationSuite, etc. when iniitalizing metadata files.
    # populate metadata files from mk_tool_directory arguments.
    raise NotImplementedError
