from pathlib import Path
from unittest import TestCase
from xd_cwl_utils.validate import validate_repo

class TestValidateRepo(TestCase):

    def test_validate_repo(self):

        print(Path.cwd())
        validate_repo()