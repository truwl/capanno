# Classes to represent metadata for command line tools.
from pathlib import Path
from abc import ABC, abstractmethod
import logging
from collections import OrderedDict
import semantic_version
from ruamel.yaml import YAML
from ruamel.yaml import safe_load
from ruamel.yaml.comments import CommentedMap
from utilities.classes.shared_properties import CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite
from utilities.get_metadata_from_biotools import make_tool_metadata_kwargs_from_biotools


object_attributes = (CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite)

class ToolMetadataBase(ABC):
    """Factor stuff out to here."""

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
        self._name = new_name


    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, new_version):
        new_version = str(new_version)
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
                raise AttributeError(f"{k} is not a valid key for {type(self)}")

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
            attr_value = getattr(self, key)
            if isinstance(attr_value, object_attributes):
                meta_map[key] = attr_value.dump()
            elif isinstance(attr_value, list):
                if isinstance(attr_value[0], object_attributes):
                    meta_map[key] = [item.dump() for item in attr_value]
                else:
                    meta_map[key] = attr_value
            else:
                meta_map[key] = attr_value
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)
        with file_path.open('w') as yaml_file:
            yaml.dump(meta_map, yaml_file)
        return


class ToolMetadata(ToolMetadataBase):
    """Class to represent metadata for a 'stand alone' command line tool."""

    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('version', '0.1.0'),  # Set to something low if not provided.
        ('identifier', None),
        ('description', None),
        ('codeRepository', CodeRepository()),
        ('license', None),
        ('WebSite', [WebSite()]),
        ('contactPoint', [Person()]),
        ('publication', [Publication()]),
        ('keywords', [Keyword()]),
        ('alternateName', None),
        ('creator', [Person()]),
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


class ParentToolMetadata(ToolMetadataBase):

    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('name', None),
        ('softwareVersion', None),
        ('version', '0.1.0'),
        ('identifier', None),
        ('featureList', None),
        ('description', None),
        ('codeRepository', CodeRepository()),
        ('license', None),
        ('WebSite', [WebSite()]),
        ('contactPoint', [Person()]),
        ('publication', [Publication()]),
        ('keywords', [Keyword()]),
        ('alternateName', None),
        ('creator', [Person()]),
        ('programmingLanguage', None),
        ('datePublished', None),
        ('downloadURL', None),
        ('extra', None)
    ])


    def make_subtool_metadata(self, subtool_name):
        if not self.featureList:
            raise ValueError(f"Cannot create subtool. featureList of {self.name} is not populated.")
        if subtool_name not in self.featureList:
            raise ValueError(f"{subtool_name} must be in the parent featureList")
        subtool_metadata = SubtoolMetadata(name=subtool_name, applicationSuite={'name': self.name, 'version': self.version, 'identifier': self.identifier})
        return subtool_metadata


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


class SubtoolMetadata(ToolMetadataBase):

    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('applicationSuite', ApplicationSuite()),
        ('name', None),
        ('version', '0.1'),
        ('description', None),
        ('keywords', Keyword()),
        ('alternateName', None),
    ])

    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict)


    @classmethod
    def initialize_from_parent(cls, parent_metadata):
        pass





