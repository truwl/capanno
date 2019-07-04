#!/usr/bin/env python3

import argparse
import sys
from utilities.validate.validate_tools import validate_parent_tool_metadata, validate_tool_metadata, validate_subtool_metadata
from utilities.validate.validate_scripts import validate_script_metadata, validate_common_script_metadata
from utilities.validate.validate_workflows import validate_workflow_metadata

parser = argparse.ArgumentParser(description="Validate metadata files.")
subparsers = parser.add_subparsers(description="Specify type of metadata to validate.", dest='command')

validate_tool = subparsers.add_parser('tool', help="Validate tool metadata.")
validate_tool.add_argument('path', help="Path to tool metadata file.")

validate_parent_tool = subparsers.add_parser('parent_tool', help="Validate parent tool metadata.")
validate_parent_tool.add_argument('path', help="Path to parent tool metadata file.")

validate_subtool = subparsers.add_parser('subtool', help="Validate subtool metadata.")
validate_subtool.add_argument('path', help="Path to subtool metadata file.")

validate_script = subparsers.add_parser('script', help="Validate script metadata")
validate_script.add_argument('path', help="Path to script metadata file.")

validate_script_common = subparsers.add_parser('script_common', help="Validate common script metadata.")
validate_script_common.add_argument('path', help="Validate common script metadata")

validate_workflow = subparsers.add_parser('workflow', help="Validate workflow metadata.")
validate_workflow.add_argument('path', help="Path to workflow metadata")

def main(args):
    if args.command == 'tool':
        validate_tool_metadata(args.path)
    elif args.command == 'parent_tool':
        validate_parent_tool_metadata(args.path)
    elif args.command == 'subtool':
        validate_subtool_metadata(args.path)
    elif args.command == 'script':
        validate_script_metadata(args.path)
    elif args.command == 'script_common':
        validate_common_script_metadata(args.path)
    elif args.command == 'workflow':
        validate_workflow_metadata(args.path)
    else:
        parser.print_help()
    return

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    sys.exit(0)