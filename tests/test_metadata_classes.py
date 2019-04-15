import os
from ruamel.yaml import safe_load
from tests.test_base import TestBase
from utilities.classes.templates import ToolMetadata, ParentToolMetadata, SubtoolMetadata


class TestMakeToolMetadata(TestBase):
    test_dict = {'name':'some_name', 'bad': 'Should not work.'}

    def test_make_tool_metadata(self):
        tm = ToolMetadata(name=TestMakeToolMetadata.test_dict['name'])
        self.assertTrue(tm.name == TestMakeToolMetadata.test_dict['name'])

    def test_bad_kwarg(self):
        with self.assertRaises(AttributeError):
            tm = ToolMetadata(bad=TestMakeToolMetadata.test_dict['bad'])

    def test_make_file(self):
        tm = ToolMetadata(name=TestMakeToolMetadata.test_dict['name'])
        test_filename = f"{TestMakeToolMetadata.test_dict['name']}--testMeta.yaml"
        tm.mk_file(test_filename)
        with open(test_filename) as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeToolMetadata.test_dict['name'])
        os.remove(test_filename)



class TestMakeSubtoolMetadata(TestBase):
    test_dict = {'name': 'subtool_name', 'bad': 'A bad key or value.'}
    def test_make_subtool_metadata(self):
        st_metadata = SubtoolMetadata(name=TestMakeSubtoolMetadata.test_dict['name'])
        self.assertTrue(st_metadata.name == TestMakeSubtoolMetadata.test_dict['name'])

    def test_make_file(self):
        st_metadata = SubtoolMetadata(name=TestMakeSubtoolMetadata.test_dict['name'])
        test_filename = f"{TestMakeSubtoolMetadata.test_dict['name']}--testMeta.yaml"
        st_metadata.mk_file(test_filename)
        with open(test_filename) as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeSubtoolMetadata.test_dict['name'])
        os.remove(test_filename)


class TestMakeParentToolMetadata(TestBase):
    test_dict = {'name': 'parent_name', 'bad': 'A bad key or value.'}

    def test_make_parent_metadata(self):
        p_metadata = ParentToolMetadata(name=TestMakeParentToolMetadata.test_dict['name'])
        self.assertTrue(p_metadata.name == TestMakeParentToolMetadata.test_dict['name'])

    def test_make_file(self):
        p_metadata = ParentToolMetadata(name=TestMakeParentToolMetadata.test_dict['name'])
        test_filename = f"{TestMakeParentToolMetadata.test_dict['name']}--testMeta.yaml"
        p_metadata.mk_file(test_filename)
        with open(test_filename) as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeParentToolMetadata.test_dict['name'])
        os.remove(test_filename)