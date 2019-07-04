
from utilities.classes.script_metadata import ScriptMetadata, CommonScriptMetadata

def validate_script_metadata(metadata_path):
    try:
        script_metadata = ScriptMetadata.load_from_file(metadata_path)
        print(f"Metadata in {metadata_path} is valid script metadata.")

    except:
        print(f"Script metadata in {metadata_path} failed validation.")
        raise
    return

def validate_common_script_metadata(metadata_path):
    try:
        common_script_metadata = CommonScriptMetadata.load_from_file(metadata_path)
        print(f"Metadata in {metadata_path} is valid common script metadata.")

    except:
        print(f"Common script metadata in {metadata_path} failed validation.")
        raise
    return