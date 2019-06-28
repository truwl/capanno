import unittest
from unittest import defaultTestLoader
from tests.test_tool_metadata_classes import TestMakeToolMetadata, TestMakeParentToolMetadata, TestMakeSubtoolMetadata
from tests.test_script_metadata import TestScriptMetadata
from tests.test_workflow_metadata_classes import TestWorkflowMetadata
from tests.test_get_metadata import TestMetadataFromBioTools


def full_suite():
    # Tool metadata classes
    suite = defaultTestLoader.loadTestsFromTestCase(TestMakeToolMetadata)
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeSubtoolMetadata))
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeParentToolMetadata))

    # Script metadata classes
    suite.addTest((defaultTestLoader.loadTestsFromTestCase(TestScriptMetadata)))

    # Workflow metadata classes
    suite.addTest((defaultTestLoader.loadTestsFromTestCase(TestWorkflowMetadata)))

    # Get metadata
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMetadataFromBioTools))
    return suite


def main():
    suite = full_suite()
    unittest.TextTestRunner().run(suite)
    return


if __name__ == '__main__':
    main()
