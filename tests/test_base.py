from unittest import TestCase
from pathlib import Path


class TestBase(TestCase):

    temp_dir = Path().cwd() / 'temp'

    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass