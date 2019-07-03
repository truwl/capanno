#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.

import os
from unittest import skip
from pathlib import Path
from ruamel.yaml import safe_load
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
        test_filename = self.temp_dir / 'script_test_metadata.yaml'
        st_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['identifier'], kwargs['identifier'])
        os.remove(test_filename)

    def test_make_script_metadata_from_file(self):
        file_path = Path().cwd() / 'cwl-scripts' / 'ENCODE_DCC' /'atac-seq-pipeline' / 'v1.1' / 'encode_ataqc-metadata.yaml'
        st_metadata = ScriptMetadata.load_from_file(file_path)
        return

    # @skip("Pass")
    def test_mk_file_with_inherited_data(self):
        file_path = file_path = Path().cwd() / 'cwl-scripts' / 'ENCODE_DCC' /'atac-seq-pipeline' / 'v1.1' / 'encode_ataqc-metadata.yaml'
        st_metadata = ScriptMetadata.load_from_file(file_path)
        test_filename = self.temp_dir / 'script2_test_metadata.yaml'
        st_metadata.mk_completed_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['license'], 'MIT')
        os.remove(test_filename)
        return

    def test_mk_file_without_inherited_data(self):
        file_path = file_path = Path().cwd() / 'cwl-scripts' / 'ENCODE_DCC' /'atac-seq-pipeline' / 'v1.1' / 'encode_ataqc-metadata.yaml'
        st_metadata = ScriptMetadata.load_from_file(file_path)
        test_filename = self.temp_dir / 'script2_test_metadata.yaml'
        st_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        with self.assertRaises(KeyError):
            license = test_file_dict['license']
        os.remove(test_filename)
        return