import os
import shutil
from unittest import skip
from pathlib import Path
from ruamel.yaml import safe_load
from tests.test_base import TestBase
from utilities.classes.tool_metadata import ToolMetadata, ParentToolMetadata, SubtoolMetadata
from utilities.add.add_tools import add_tool, add_subtool, add_parent_tool


class TestMakeToolMetadata(TestBase):
    test_dict = {'name':'some_name', 'bad': 'Should not work.', 'softwareVersion': 1}

    def test_make_tool_metadata(self):
        tm = ToolMetadata(name=TestMakeToolMetadata.test_dict['name'], softwareVersion=TestMakeToolMetadata.test_dict['softwareVersion'])
        self.assertTrue(tm.name == TestMakeToolMetadata.test_dict['name'])

    def test_bad_kwarg(self):
        with self.assertRaises(AttributeError):
            tm = ToolMetadata(bad=TestMakeToolMetadata.test_dict['bad'])

    def test_make_file(self):
        tm = ToolMetadata(name=TestMakeToolMetadata.test_dict['name'], softwareVersion=0.1)
        test_filename = self.temp_dir() / 'tool_test.yaml'
        tm.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeToolMetadata.test_dict['name'])
        os.remove(test_filename)

    def test_make_from_biotools(self):
        biotools_id = 'star'
        file_name = self.temp_dir() / 'biotools_test.yaml'
        biotools_meta = ToolMetadata.create_from_biotools(biotools_id, softwareVersion=1)
        biotools_meta.mk_file(file_name)
        os.remove(file_name)
        return


@skip('Cant make SubtoolMetadata from kwargs for now. Needs parentMetadata')
class TestMakeSubtoolMetadata(TestBase):
    test_dict = {'name': 'subtool_name', 'bad': 'A bad key or value.'}
    def test_make_subtool_metadata(self):
        st_metadata = SubtoolMetadata(name=TestMakeSubtoolMetadata.test_dict['name'])
        self.assertTrue(st_metadata.name == TestMakeSubtoolMetadata.test_dict['name'])

    def test_make_file(self):
        st_metadata = SubtoolMetadata(name=TestMakeSubtoolMetadata.test_dict['name'])
        test_filename = self.temp_dir() / 'subtool_test.yaml'
        st_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeSubtoolMetadata.test_dict['name'])
        os.remove(test_filename)


# @skip('Pass to isolate tests')
class TestMakeParentToolMetadata(TestBase):
    test_dict = {'name': 'parent_name', 'bad': 'A bad key or value.'}

    def test_make_parent_metadata(self):
        p_metadata = ParentToolMetadata(name=TestMakeParentToolMetadata.test_dict['name'], softwareVersion=1)
        self.assertTrue(p_metadata.name == TestMakeParentToolMetadata.test_dict['name'])

    def test_make_file(self):
        p_metadata = ParentToolMetadata(name=TestMakeParentToolMetadata.test_dict['name'], softwareVersion=1)

        test_filename = self.temp_dir() / 'parent.yaml'
        p_metadata.mk_file(test_filename)
        with test_filename.open('r') as file:
            test_file_dict = safe_load(file)
        self.assertEqual(test_file_dict['name'], TestMakeParentToolMetadata.test_dict['name'])
        os.remove(test_filename)

class TestAddTools(TestBase):

    _foo_tool_name = "ReallyRidiculousToolName"  # Make variable name different so don't overwrite or inadvertently access.
    _foo_softwareVersion = "ThisIsAWeirdSoftwareVersion"
    _foo_subtool_name = 'YeahNo'

    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        test_dir = Path().cwd() / 'cwl-tools' / self._foo_tool_name
        shutil.rmtree(test_dir)  # dangerous. Make sure cls.tool_name isn't messed with.


    def test_add_tool(self):
        new_tool_path = add_tool(self._foo_tool_name, self._foo_softwareVersion)
        return

    def test_add_parent_tool(self):
        new_tool_path = add_parent_tool(self._foo_tool_name, self._foo_softwareVersion, self._foo_subtool_name)
        return

    def test_add_subtool(self):
        parent_tool_path = add_parent_tool(self._foo_tool_name, self._foo_softwareVersion, self._foo_subtool_name)
        subtool_path = add_subtool(parent_tool_path, self._foo_subtool_name)
        return

