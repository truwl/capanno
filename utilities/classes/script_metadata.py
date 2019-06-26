#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.



from pathlib import Path
from collections import OrderedDict
from ruamel.yaml import YAML, safe_load
from ruamel.yaml.comments import CommentedMap
from utilities.classes.shared_properties import WebSite, CodeRepository, Person, Publication, Keyword
from utilities.classes.metadata_base import MetadataBase

class ScriptMetadataBase(MetadataBase):
    """Factor stuff out to here."""

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        # Can put validators here.
        self._name = new_name


class ScriptMetadata(ScriptMetadataBase):

    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('description', None),
        ('identifier', None),
        ('version', '0.1'),
        ('WebSite', [WebSite()]),
        ('codeRepository', CodeRepository()),
        ('license', None),
        ('contactPoint', [Person()]),
        ('publication', [Publication()]),
        ('keywords', [Keyword()]),
        ('alternateName', None),
        ('creator', [Person()]),
        ('programmingLanguage', None),
        ('datePublished', None),
        ('downloadURL', None),
        ('parentScripts', {'name': None, 'version': None, 'identifier': None}),
        ('tools', None),
        ('commonMetadataPath', None),
        ('parentMetadata', None),
    ])

    def _load_common_metadata(self, file_path):
        # Start with returning a list of dicts.
        if self.commonMetadataPath:
            common_meta_list = []
            dirname = file_path.parents[0]
            for rel_path in self.commonMetadataPath:
                full_path = dirname / rel_path
                with full_path.resolve().open('r') as f:
                    common_meta_dict = safe_load(f)
                common_meta_list.append(ScriptMetadata(**common_meta_dict))
            return common_meta_list
        return



    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        new_instance = cls(**file_dict)
        new_instance.parentMetadata = new_instance._load_common_metadata(file_path)
        return new_instance

    def dump_with_parent_data(self):
        raise NotImplementedError
