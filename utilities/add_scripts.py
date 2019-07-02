
from pathlib import Path
from utilities.classes.script_metadata import ScriptMetadata, CommonScriptMetadata

def _get_script_directory(group_name, project_name, software_version):
    base_path = Path().cwd() / 'cwl-scripts' / group_name / project_name / software_version
    if base_path.exists():
        if not base_path.is_dir():
            raise TypeError(f"{base_path} must be a directory")
    else:
        base_path.mkdir(parents=True)
    return base_path


def add_common_script_metadata(group_name, project_name, software_version, filename, **kwargs):
    path = _get_script_directory(group_name, project_name, software_version) / 'common'
    path.mkdir(exist_ok=True)
    file_path = path / filename
    script_metadata = CommonScriptMetadata(**kwargs)
    script_metadata.mk_file(file_path)
    return


def add_script(script_name, script_version, parent_metadata=None):
    script_version = str(script_version)
