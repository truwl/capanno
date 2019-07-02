#!/usr/bin/env python3

# Module to call from command line.

import sys
import argparse
from utilities.add_tools import add_parent_tool, add_tool, add_subtool

parser = argparse.ArgumentParser(description='Initialize metadata and directories for a tool, script, or workflow.')
subparsers = parser.add_subparsers(description='Specify the command to run.', dest='command')

# add_tool parser
addtool = subparsers.add_parser('add_tool', help='add a new standalone tool.')
addtool.add_argument('tool_name', type=str, help="The name of the tool to add.")
addtool.add_argument('tool_version', type=str, help="The version of the tool to add.")
addtool.add_argument('--biotoolsID', type=str, help='biotools id from https://bio.tools')

# add_parent_tool parser
add_parent = subparsers.add_parser('add_parent_tool', help='add parent tool/common metadata')
add_parent.add_argument('tool_name', type=str, help='The name of the tool to add.')
add_parent.add_argument('tool_version', type=str, help="The version of the tool to add.")
add_parent.add_argument('subtools', help='list of subtools of the parent tool.', nargs='+')
add_parent.add_argument('--biotoolsID', type=str, help='biotools id from bio.tools')

# add_subtool parser
addsubtool = subparsers.add_parser('add_subtool', help='add a subtool. A parent must exist.')
addsubtool.add_argument('parent_path', help='The absolute or relative path of the metadata file that describes the parent of the subtool.')
addsubtool.add_argument('subtool_name', help="The name of the subtool. The subtool name must be present in the parent's featureList field.")




def main(args):
    if args.command == 'add_tool':
        add_tool(args.tool_name, args.tool_version, args.biotoolsID)
    elif args.command == 'add_parent_tool':
        add_parent_tool(args.tool_name, args.tool_version, args.subtools, args.biotoolsID)
    elif args.command == 'add_subtool':
        add_subtool(args.parent_path, args.subtool_name)
    else:
        raise NotImplementedError(f"Deal with command {args.command}")
    return

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    sys.exit(0)