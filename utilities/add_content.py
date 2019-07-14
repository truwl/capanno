#!/usr/bin/env python3

# Module to call from command line.

import sys
import argparse
from utilities.add.add_tools import add_parent_tool, add_tool, add_subtool
from utilities.add.add_scripts import add_script, add_common_script_metadata
from utilities.add.add_workflows import add_workflow

parser = argparse.ArgumentParser(description='Initialize metadata and directories for a tool, script, or workflow.')
subparsers = parser.add_subparsers(description='Specify the command to run.', dest='command')

# add_tool parser
addtool = subparsers.add_parser('tool', help='add a new standalone tool.')
addtool.add_argument('tool_name', type=str, help="The name of the tool to add.")
addtool.add_argument('tool_version', type=str, help="The version of the tool to add.")
addtool.add_argument('--biotoolsID', type=str, help='biotools id from https://bio.tools')

# add_parent_tool parser
add_parent = subparsers.add_parser('parent_tool', help='add parent tool/common metadata')
add_parent.add_argument('tool_name', type=str, help='The name of the tool to add.')
add_parent.add_argument('tool_version', type=str, help="The version of the tool to add.")
add_parent.add_argument('subtools', help='list of subtools of the parent tool.', nargs='+')
add_parent.add_argument('--biotoolsID', type=str, help='biotools id from bio.tools')

# add_subtool parser
addsubtool = subparsers.add_parser('subtool', help='add a subtool. A parent must exist.')
addsubtool.add_argument('parent_path', help='The absolute or relative path of the metadata file that describes the parent of the subtool.')
addsubtool.add_argument('subtool_name', help="The name of the subtool. The subtool name must be present in the parent's featureList field.")

# add_common_script_parser
addscriptcommon = subparsers.add_parser('common_script', help='add script metadata that other scripts can inherit from')
addscriptcommon.add_argument('group_name', help='The name of the group directory that the script will go into.')
addscriptcommon.add_argument('project_name', help='The name of the project directory that the script will go into.')
addscriptcommon.add_argument('script_version', help='The version of the script.')
addscriptcommon.add_argument('filename', help="The name of common metadata file. File extensions and '-metadata.yaml' postfix must be omitted.")

# add_script parser
addscript = subparsers.add_parser('script', help='add a script')
addscript.add_argument('group_name', help='The name of the group directory that the script will go into.')
addscript.add_argument('project_name', help='The name of the project directory that the script will go into.')
addscript.add_argument('script_version', help='The version of the script.')
addscript.add_argument('script_name', help='The name of the script. File extensions should be omitted. Replace spaces with underscores')
addscript.add_argument('--parent_metadata', '-p', nargs='+', help="path(s) to common script metadata that the script should inherit from.")
# Could add flags for more kwargs here. Probably easiest to just add to file though.

addworkflow = subparsers.add_parser('workflow', help='add a workflow')
addworkflow.add_argument('group_name', help='The name of the group directory that the workflow will go into.')
addworkflow.add_argument('workflow_name', help='The name of the workflow. Replace spaces with underscores.')
addworkflow.add_argument('workflow_version', help='The version of the workflow.')


def main(args):
    if args.command == 'tool':
        add_tool(args.tool_name, args.tool_version, args.biotoolsID)
    elif args.command == 'parent_tool':
        add_parent_tool(args.tool_name, args.tool_version, args.subtools, args.biotoolsID)
    elif args.command == 'subtool':
        add_subtool(args.parent_path, args.subtool_name)
    elif args.command == 'common_script':
        add_common_script_metadata(args.group_name, args.project_name, args.script_version, args.filename)
    elif args.command == 'script':
        add_script(args.group_name, args.project_name, args.script_version, args.script_name, parentMetadata=args.parent_metadata)
    elif args.command == 'workflow':
        add_workflow(args.group_name, args.workflow_name, args.workflow_version)
    else:
        parser.print_help()
    return

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    sys.exit(0)