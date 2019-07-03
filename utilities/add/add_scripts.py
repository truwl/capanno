
from pathlib import Path
from utilities.classes.script_metadata import ScriptMetadata, CommonScriptMetadata

def _get_script_directory(group_name, project_name, script_version):
    """
    Gets the script directory specified. Makes the directory if it doesn't exist.
    :param group_name:
    :param project_name:
    :param script_version:
    :return:
    """
    base_path = Path().cwd() / 'cwl-scripts' / group_name / project_name / script_version
    if base_path.exists():
        if not base_path.is_dir():
            raise TypeError(f"{base_path} must be a directory")
    else:
        base_path.mkdir(parents=True)
    return base_path


def add_common_script_metadata(group_name, project_name, script_version, filename, **kwargs):
    path = _get_script_directory(group_name, project_name, script_version) / 'common'
    path.mkdir(exist_ok=True)
    filename = f"{filename}-metadata.yaml"
    file_path = path / filename
    script_metadata = CommonScriptMetadata(**kwargs)
    script_metadata.mk_file(file_path)
    return


def add_script(group_name, project_name, script_version, script_name, **kwargs):
    script_version = str(script_version)
    script_dir = _get_script_directory(group_name, project_name, script_version) / script_name
    script_dir.mkdir()
    script_metadata = ScriptMetadata(name=script_name, softwareVersion=script_version, **kwargs)
    filename = f"{script_name}-metadata.yaml"
    script_metadata.mk_file(script_dir / filename)
    instances_dir = script_dir / 'instances'
    instances_dir.mkdir()
    return


