
from tests.test_base import TestBase
from utilities.content_maps import make_tools_map, make_script_map, make_script_maps

class TestToolMaps(TestBase):

    def test_make_tools_map(self):
        make_tools_map()
        return


    def test_make_script_map(self):
        group_name = 'ENCODE-DCC'
        project = 'atac-seq-pipeline'
        version = 'v1.1'
        make_script_map(group_name, project, version)

    def test_make_script_maps(self):
        make_script_maps()
        return