#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.

import os
from pathlib import Path
from ruamel.yaml import safe_load
from tests.test_base import TestBase
from utilities.classes.workflow_metadata import WorkflowMetadata


class TestWorkflowMetadata(TestBase):

    def test_make_workflow_metadata(self):
        test_name = 'Test wf name'
        wf = WorkflowMetadata(name=test_name)
        self.assertEqual(test_name, wf.name)

    def test_load_from_file(self):
        file_path = Path().cwd() / 'cwl-workflows' / 'example_workflows' / 'cat_sort' / '1.0' / 'cat_sort-metadata.yaml'
        wf_meta = WorkflowMetadata.load_from_file(file_path)
        self.assertEqual(wf_meta.name, 'cat sort')

    def test_mk_file(self):
        file_path = Path().cwd() / 'cwl-workflows' / 'example_workflows' / 'cat_sort' / '1.0' / 'cat_sort-metadata.yaml'
        wf_meta = WorkflowMetadata.load_from_file(file_path)
        test_filename = self.temp_dir / 'wf_test_metadata.yaml'
        wf_meta.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['identifier'], 'WF_1f4d8f.cb')
        os.remove(test_filename)