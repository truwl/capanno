from pathlib import Path
from unittest import TestCase
from xd_cwl_utils.validate import validate_repo

class TestValidateRepo(TestCase):

    def test_validate_repo(self):

        base_dir = Path.cwd()
        validate_repo(base_dir=base_dir)