
from pathlib import Path
from tests.test_base import TestBase
from xd_cwl_utils.content_maps import make_tools_map, make_script_map, make_script_maps
class TestToolMaps(TestBase):

    def test_make_tools_map(self):
        tool_dir = Path.cwd() / 'cwl-tools'
        outfile = Path.cwd() / 'content_maps' / 'tool-marps.yaml'
        make_tools_map(tool_dir, outfile)
        return


    def test_make_script_map(self):
        group_name = 'ENCODE-DCC'
        project = 'atac-seq-pipeline'
        version = 'v1.1'
        make_script_map(group_name, project, version)

    def test_make_script_maps(self):
        make_script_maps()
        return