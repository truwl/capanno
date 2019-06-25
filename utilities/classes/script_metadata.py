#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.


from abc import ABC, abstractclassmethod
from collections import OrderedDict
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

from utilities.classes.tool_metadata import ToolMetadataBase

class ScriptMetadataBase(ABC):
    pass


class ScriptMetadata(ScriptMetadataBase):
    _init_metadata = OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('description', None),
        ('identifier', None),
        ('version', None),
        ('WebSite', [{'name': None, 'description': None, 'URL': None}]),
        ('codeRepository', dict([('name', None), ('URL', None)])),
        ('license', None),
        ('contactPoint', [{'name': None, 'email': None, 'identifier': None}]),
        ('publication', [{'identifier': None, 'headline': None}]),
        ('keywords', None),
        ('alternateName', None),
        ('creator', [{'name': None, 'email': None, 'identifier': None}]),
        ('programmingLanguage', None),
        ('datePublished', None),
        ('downloadURL', None),
        ('parentScripts', {'name': None, 'version': None, 'identifier': None}),
        ('parentMetadata', None),
    ])

    @classmethod
    def _metafile_keys(cls):
        return list(cls._init_metadata.keys())

    def __init__(self, **kwargs):
        super().__init__()
        for k, v in self._init_metadata.items():
            setattr(self, k, v)

        for k, v in kwargs.items():
            if not k in self._init_metadata:
                raise KeyError(f"{k} is not a valid key for ScriptMetadata")
            setattr(self, k, v)
        return

    def mk_file(self, file_path):
        keys = ScriptMetadata._metafile_keys()
        meta_map = CommentedMap()
        for key in keys:
            meta_map[key] = getattr(self, key)
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with open(file_path, 'w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return