#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.

import os
from unittest import skip
from pathlib import Path
from ruamel.yaml import safe_load
from content_maps import script_maps
from tests.test_base import TestBase
from utilities.classes.script_metadata import ScriptMetadata

class TestScriptMetadata(TestBase):

    def test_make_script_metadata_from_kwargs(self):
        kwargs = {'name': 'test_script', 'softwareVersion': 1, 'identifier': 'ST_abcdef.12'}
        st_metadata = ScriptMetadata(**kwargs)
        self.assertEqual(st_metadata.name, kwargs['name'])
        self.assertEqual(st_metadata.identifier, kwargs['identifier'])
        return

    def test_mk_file_from_script_metadata(self):
        kwargs = {'name': 'test_script', 'softwareVersion': 1, 'identifier': 'ST_abcdef.12'}
        st_metadata = ScriptMetadata(**kwargs)
        test_filename = Path(os.environ.get('TEST_TMP_DIR')) / 'script_test_metadata.yaml'
        st_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['identifier'], kwargs['identifier'])
        os.remove(test_filename)

    def test_make_script_metadata_from_file(self):
        metadata_path = TestBase.get_metadata_path(script_maps.ENCODE_atac_seq['ST_43baaf.f7']) # encode_ataqc.cwl
        st_metadata = ScriptMetadata.load_from_file(metadata_path)
        return

    # @skip("Pass")
    def test_mk_file_with_inherited_data(self):
        metadata_path = TestBase.get_metadata_path(script_maps.ENCODE_atac_seq['ST_43baaf.f7'])
        st_metadata = ScriptMetadata.load_from_file(metadata_path)
        test_filename = Path(os.environ.get('TEST_TMP_DIR')) / 'script2_test_metadata.yaml'
        st_metadata.mk_completed_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['license'], 'MIT')
        os.remove(test_filename)
        return

    def test_mk_file_without_inherited_data(self):
        metadata_path = TestBase.get_metadata_path(script_maps.ENCODE_atac_seq['ST_43baaf.f7'])
        st_metadata = ScriptMetadata.load_from_file(metadata_path)
        test_filename = Path(os.environ.get('TEST_TMP_DIR')) / 'script3_test_metadata.yaml'
        st_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        with self.assertRaises(KeyError):
            license = test_file_dict['license']
        os.remove(test_filename)
        return