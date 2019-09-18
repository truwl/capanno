import unittest
import logging
import tempfile
import os
from unittest import defaultTestLoader
from xd_cwl_utils.validate import validate_repo

class TestValidateRepo(unittest.TestCase):

    def test_validate_repo(self):
        validate_repo()


def test_suite():
    suite = defaultTestLoader.loadTestsFromTestCase(TestValidateRepo)
    return suite


def main():
    logging.basicConfig(level=logging.CRITICAL)
    os.environ['CONFIG_KEY'] = 'DEFAULT'
    suite = test_suite()
    unittest.TextTestRunner().run(suite)
    return


if __name__ == '__main__':
    main()
