
import logging
import os
from unittest import defaultTestLoader, TextTestRunner
from tests.test_validate_repo import TestValidateRepo


def test_suite():
    suite = defaultTestLoader.loadTestsFromTestCase(TestValidateRepo)
    return suite


def main():
    logging.basicConfig(level=logging.CRITICAL)
    os.environ['CONFIG_KEY'] = 'DEFAULT'
    suite = test_suite()
    TextTestRunner().run(suite)
    return


if __name__ == '__main__':
    main()
