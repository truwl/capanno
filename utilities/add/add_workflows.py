
from pathlib import Path
from utilities.classes.workflow_metadata import WorkflowMetadata


def _get_workflow_directory(group_name, workflow_name, workflow_version):
    """
    Gets the script directory specified. Makes the directory if it doesn't exist.
    :param group_name:
    :param project_name:
    :param script_version:
    :return:
    """
    base_path = Path().cwd() / 'cwl-workflows' / group_name / workflow_name / workflow_version
    if base_path.exists():
        if not base_path.is_dir():
            raise TypeError(f"{base_path} must be a directory")
    else:
        base_path.mkdir(parents=True)
    return base_path


def add_workflow(group_name, workflow_name, workflow_version, **kwargs ):
    workflow_version = str(workflow_version)
    wf_path = _get_workflow_directory(group_name, workflow_name, workflow_version)
    wf_metadata = WorkflowMetadata(name=workflow_name, softwareVersion=workflow_version, **kwargs)
    filename = f"{workflow_name}-metadata.yaml"
    wf_metadata.mk_file(wf_path / filename)
    instances_dir = wf_path / 'instances'
    instances_dir.mkdir()
    return