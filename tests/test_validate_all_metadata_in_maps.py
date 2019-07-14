from tests.test_base import TestBase
from pathlib import Path
from ruamel.yaml import safe_load
from content_maps import tool_maps, script_maps, workflow_maps
from utilities.validate import main

class TestValidateContent(TestBase):

    def test_validate_tools(self):
        for identifier, tool_path in tool_maps.everything.items():
            path = Path(tool_path)
            if path.suffix == '.yaml':
                if not 'common' in path.parts:
                    raise ValueError(f"Have a parent-like path that is not in a common directory {path}")
                meta_type = 'parent_tool'
                meta_path = path
            else:  # either a subtool or standalone tool.
                meta_path = TestBase.get_metadata_path(path)
                with meta_path.open('r') as f:
                    meta_dict = safe_load(f)
                if meta_dict.get('parentMetadata'):
                    meta_type = 'subtool'
                else:
                    meta_type = 'tool'
            main(meta_type, meta_path)
        return

    def test_validate_scripts(self):
        for script_identifier, script_path in script_maps.everything.items():
            metadata_path = TestBase.get_metadata_path(script_path)
            main('script', metadata_path)
        return


    def test_validate_workflows(self):
        for workflow_identifier, workflow_cwl_path in workflow_maps.everything.items():
            metadata_path = TestBase.get_metadata_path(workflow_cwl_path)
            main('workflow', metadata_path)
        return