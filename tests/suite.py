import unittest
import logging
from unittest import defaultTestLoader
from tests.test_tool_metadata import TestMakeToolMetadata, TestMakeParentToolMetadata, TestMakeSubtoolMetadata, TestAddTools
from tests.test_script_metadata import TestScriptMetadata
from tests.test_workflow_metadata import TestWorkflowMetadata
from tests.test_get_metadata import TestMetadataFromBioTools
from tests.test_validate import TestValidateMetadata




def full_suite():
    # Tool metadata classes
    suite = defaultTestLoader.loadTestsFromTestCase(TestMakeToolMetadata)
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeSubtoolMetadata))
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeParentToolMetadata))
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestAddTools))

    # Script metadata classes
    suite.addTest((defaultTestLoader.loadTestsFromTestCase(TestScriptMetadata)))

    # Workflow metadata classes
    suite.addTest((defaultTestLoader.loadTestsFromTestCase(TestWorkflowMetadata)))

    # Get metadata
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMetadataFromBioTools))

    # Validate
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestValidateMetadata))

    return suite



def main():
    logging.basicConfig(level=logging.CRITICAL)
    suite = full_suite()
    unittest.TextTestRunner().run(suite)
    return


if __name__ == '__main__':
    main()
