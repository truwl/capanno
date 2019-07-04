#!/usr/bin/env python3

import argparse
import sys
from utilities.validate.validate_tools import validate_parent_tool_metadata, validate_tool_metadata, validate_subtool_metadata

parser = argparse.ArgumentParser(description="Validate metadata files.")
subparsers = parser.add_subparsers(description="Specify type of metadata to validate.", dest='command')

validate_tool = subparsers.add_parser('tool', help="Validate tool metadata.")
validate_tool.add_argument('path', help="Path to tool metadata file.")

validate_parent_tool = subparsers.add_parser('parent_tool', help="Validate parent tool metadata.")
validate_parent_tool.add_argument('path', help="Path to parent tool metadata file.")

validate_subtool = subparsers.add_parser('subtool', help="Validate subtool metadata.")
validate_subtool.add_argument('path', help="Path to subtool metadata file.")

def main(args):
    if args.command == 'tool':
        validate_tool_metadata(args.path)
    elif args.command == 'parent_tool':
        validate_parent_tool_metadata(args.path)
    elif args.command == 'subtool':
        validate_subtool_metadata(args.path)
    else:
        parser.print_help()
    return

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    sys.exit(0)