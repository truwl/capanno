# Classes to represent metadata for command line tools.
from pathlib import Path
from abc import ABC, abstractmethod
import logging
from collections import OrderedDict
import semantic_version
from ruamel.yaml import YAML
from ruamel.yaml import safe_load
from ruamel.yaml.comments import CommentedMap
from utilities.get_metadata_from_biotools import make_tool_metadata_kwargs_from_biotools


class MetadataBase(ABC):
    """Will factor stuff out to here eventually."""


    @staticmethod
    @abstractmethod
    def _init_metadata():
        return OrderedDict([('name', None), ('version', None)])


    # Properties/Fields that every metadata class will have

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        # Can put validators here.
        print('hit it.')
        self._name = new_name


    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, new_version):
        is_semantic = semantic_version.validate(new_version)
        if is_semantic:
            v = semantic_version.Version(new_version)
        else:
            logging.info(f"'{new_version}' is not a properly formatted semantic version")
            try:
                v = semantic_version.Version(new_version, partial=True)
            except ValueError:
                raise ValueError(f"'{new_version}'' is not a valid partial semantic version.")
            v = semantic_version.Version.coerce(str(v))
        self._version = str(v)



    def _get_metafile_keys(self):
        return list(self._init_metadata())


    def __init__(self, **kwargs):


        metadata = self._init_metadata()

        for k, v in kwargs.items():
            if not k in metadata:
                raise AttributeError(f"{k} is not a valid key for ToolMetadata")

        for k, v in metadata.items():
            if k in kwargs:
                setattr(self, k, kwargs[k])
            else:
                setattr(self, k, v)
        return

    def mk_file(self, file_path):
        file_path = Path(file_path)
        meta_map = CommentedMap()
        keys = self._get_metafile_keys()
        for key in keys:
            meta_map[key] = getattr(self, key)
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with file_path.open('w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return



class CodeRepository:
    pass

class WebSite:
    pass

class Publication:
    pass

class Person:
    pass



class ToolMetadata(MetadataBase):
    """Class to represent metadata for a command line tool."""


    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('version', None),
        ('description', None),
        ('codeRepository', dict([('name', None), ('URL', None)])),
        ('license', None),
        ('WebSite', [{'name': None, 'description': None, 'URL': None}]),
        ('contactPoint', [{'name': None, 'email': None, 'identifier': None}]),
        ('publication', [{'identifier': None, 'headline': None}]),
        ('keywords', None),
        ('alternateName', None),
        ('creator', [{'name': None, 'email': None, 'identifier': None}]),
        ('programmingLanguage', None),
        ('datePublished', None),
        ('downloadURL', None),
        ('extra', None),
    ])


    # Class factory methods

    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict)

    @classmethod
    def create_from_biotools(cls, biotools_id, version='0.1.1'):
        kwargs = make_tool_metadata_kwargs_from_biotools(biotools_id)
        return cls(**kwargs, version=version)



class SubtoolMetadata(MetadataBase):
    _init_metadata = OrderedDict([
        ('applicationSuite', {'name': None, 'softwareVersion': None, 'identifier': None}),
        ('name', None),
        ('version', None),
        ('description', None),
        ('keywords', None),
        ('alternateName', None),
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
                raise KeyError(f"{k} is not a valid key for ToolMetadata")
            setattr(self, k, v)
        return

    @classmethod
    def initialize_from_parent(cls, parent_metadata):
        pass

    def mk_file(self, file_path):
        keys = SubtoolMetadata._metafile_keys()
        meta_map = CommentedMap()
        for key in keys:
            meta_map[key] = getattr(self, key)
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with open(file_path, 'w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return


class ParentToolMetadata(MetadataBase):
    _init_metadata = OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('featureList', None),
        ('description', None),
        ('codeRepository', dict([('name', None), ('URL', None)])),
        ('license', None),
        ('WebSite', [{'name': None, 'description': None, 'URL': None}]),
        ('contactPoint', [{'name': None, 'email': None, 'identifier': None}]),
        ('publication', [{'identifier': None, 'headline': None}]),
        ('keywords', None),
        ('alternateName', None),
        ('creator', [{'name': None, 'email': None, 'identifier': None}]),
        ('programmingLanguage', None),
        ('datePublished', None),
        ('downloadURL', None)
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
                raise KeyError(f"{k} is not a valid key for ToolMetadata")
            setattr(self, k, v)
        return

    def mk_file(self, file_path):
        keys = ParentToolMetadata._metafile_keys()
        meta_map = CommentedMap()
        for key in keys:
            meta_map[key] = getattr(self, key)
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with open(file_path, 'w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return


class ScriptMetadata(MetadataBase):
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

class WorkflowMetadata(MetadataBase):
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
        ('callMap', [{'id': None, 'identifier': None}])
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
                raise KeyError(f"{k} is not a valid key for WorkflowMetadata")
            setattr(self, k, v)
        return

    def mk_file(self, file_path):
        keys = WorkflowMetadata._metafile_keys()
        meta_map = CommentedMap()
        for key in keys:
            meta_map[key] = getattr(self, key)
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with open(file_path, 'w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return
