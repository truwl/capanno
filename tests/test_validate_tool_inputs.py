

from tests.test_base import TestBase
from utilities.validate_tool_inputs import validate_tool_inputs

class TestValidateInputs(TestBase):

    def test_validate_tool_inputs(self):
        tool_name = 'samtools'
        tool_version = '1.3'
        subtool = 'flagstat'
        instance_hash = '395d'
        validate_tool_inputs(tool_name, tool_version, instance_hash, subtool)
        return
