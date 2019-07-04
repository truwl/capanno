
from utilities.classes.tool_metadata import ToolMetadata, ParentToolMetadata, SubtoolMetadata

def validate_tool_metadata(metadata_path):
    try:
        tool_metadata = ToolMetadata.load_from_file(metadata_path)
        print(f"Metadata im {metadata_path} is valid tool metadata.")
    except:
        print(f"Tool metadata in {metadata_path} failed validation")
        raise
    return


def validate_parent_tool_metadata(metadata_path):
    try:
        parent_metadata = ParentToolMetadata.load_from_file(metadata_path)
        print(f"Metadata im {metadata_path} is valid parent tool metadata.")
    except:
        print(f"Parent tool metadata in {metadata_path} failed validation")
        raise
    return

def validate_subtool_metadata(metadata_path):
    try:
        subtool_metadata = SubtoolMetadata.load_from_file(metadata_path)
        print(f"Metadata im {metadata_path} is valid subtool metadata.")
    except:
        print(f"Subtool metadata in {metadata_path} failed validation")
        raise
    return