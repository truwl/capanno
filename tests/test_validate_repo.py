from unittest import TestCase
from xd_cwl_utils.validate import validate_repo

class TestValidateRepo(TestCase):

    def test_validate_repo(self):
        validate_repo()