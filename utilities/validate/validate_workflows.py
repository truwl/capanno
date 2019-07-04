
from utilities.classes.workflow_metadata import WorkflowMetadata

def validate_workflow_metadata(metadata_path):
    try:
        wf_metadata = WorkflowMetadata.load_from_file(metadata_path)
        print(f"Metadata in {metadata_path} is valid workflow metadata.")
    except:
        print(f"Workflow metadata in {metadata_path} failed validation.")
        raise
    return