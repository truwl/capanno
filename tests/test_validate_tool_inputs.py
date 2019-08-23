

from tests.test_base import TestBase
from utilities.validate_tool_inputs import validate_tool_inputs

class TestValidateInputs(TestBase):

    def test_validate_tool_inputs(self):
        tool_name = 'cat'
        tool_version = '8.25'
        subtool = None
        instance_hash = '8a6c'
        validate_tool_inputs(tool_name, tool_version, instance_hash, subtool)
        return

