import os
from tests.test_base import TestBase
from ruamel.yaml import safe_load


class TestMetadataFromBioTools(TestBase):

    class_data = {'biotoolsID': 'signalP',
                  'name': 'SignalP'}

    def test_populate_file(self):
        pass
        # biotools_ID = TestMetadataFromBioTools.class_data['biotoolsID']
        # file_name = make_metadata_file_from_biotools(biotools_ID)
        # with open(file_name, 'r') as file:
        #     meta_file_dict = safe_load(file)
        # self.assertEqual(meta_file_dict['name'], TestMetadataFromBioTools.class_data['name'])
        # os.remove(file_name)
