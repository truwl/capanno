from unittest import TestCase
from pathlib import Path


class TestBase(TestCase):

    _temp_dir = Path().cwd() / 'temp'

    @classmethod
    def temp_dir(cls):
        return cls._temp_dir

    @staticmethod
    def get_metadata_path(rel_primary_file_path):
        full_path = Path.cwd() / rel_primary_file_path
        path_dir = full_path.parent
        file_name = full_path.stem
        metafile_name = file_name + '-metadata.yaml'
        metafile_path = path_dir / metafile_name
        return metafile_path



    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass