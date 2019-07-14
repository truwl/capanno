#!/usr/bin/env python3

import argparse
import sys
import logging
from utilities.classes.tool_metadata import ToolMetadata, ParentToolMetadata, SubtoolMetadata
from utilities.classes.script_metadata import ScriptMetadata, CommonScriptMetadata
from utilities.classes.workflow_metadata import WorkflowMetadata


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

def metadata_validator_factory(class_to_validate):
    def metadata_validator(metadata_path):
        try:
            metadata = class_to_validate.load_from_file(metadata_path)
            # print(f"Metadata in {metadata_path} is valid {str(class_to_validate)}")
            logging.info(f"Metadata in {metadata_path} is valid {str(class_to_validate)}")
        except:
            logging.error(f"{str(class_to_validate)} in {metadata_path} failed validation.")
            raise
        return
    return metadata_validator

def main(command, path):
    if command == 'tool':
        validate_tool_metadata = metadata_validator_factory(ToolMetadata)
        validate_tool_metadata(path)
    elif command == 'parent_tool':
        validate_parent_tool_metadata = metadata_validator_factory(ParentToolMetadata)
        validate_parent_tool_metadata(path)
    elif command == 'subtool':
        validate_subtool_metadata = metadata_validator_factory(SubtoolMetadata)
        validate_subtool_metadata(path)
    elif command == 'script':
        validate_script_metadata = metadata_validator_factory(ScriptMetadata)
        try:
            validate_script_metadata(path)
        except:
            logging.error(f"Problem with {path}")
            raise
    elif command == 'script_common':
        validate_common_script_metadata = metadata_validator_factory(CommonScriptMetadata)
        validate_common_script_metadata(path)
    elif command == 'workflow':
        validate_workflow_metadata = metadata_validator_factory(WorkflowMetadata)
        validate_workflow_metadata(path)
    else:
        parser.print_help()
    return

if __name__ == "__main__":
    args = parser.parse_args()
    main(args.command, args.path)
    sys.exit(0)