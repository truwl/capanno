import unittest
from unittest import defaultTestLoader
from tests.test_tool_metadata_classes import TestMakeToolMetadata, TestMakeParentToolMetadata, TestMakeSubtoolMetadata
from tests.test_get_metadata import TestMetadataFromBioTools


def full_suite():
    # Metadata classes
    suite = defaultTestLoader.loadTestsFromTestCase(TestMakeToolMetadata)
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeSubtoolMetadata))
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMakeParentToolMetadata))

    # Get metadata
    suite.addTest(defaultTestLoader.loadTestsFromTestCase(TestMetadataFromBioTools))
    return suite


def main():
    suite = full_suite()
    unittest.TextTestRunner().run(suite)
    return


if __name__ == '__main__':
    main()
