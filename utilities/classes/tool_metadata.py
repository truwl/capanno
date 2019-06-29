# Classes to represent metadata for command line tools.
from pathlib import Path
from collections import OrderedDict
from ruamel.yaml import safe_load
from utilities.classes.metadata_base import MetadataBase
from utilities.classes.shared_properties import CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite, \
    IOObjectItem
from utilities.get_metadata_from_biotools import make_tool_metadata_kwargs_from_biotools



class ToolMetadataBase(MetadataBase):
    """Factor stuff out to here."""


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        # Can put validators here.
        self._name = new_name





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
    def initialize_from_parent(cls, parent_metadata, subtool_name):
        subtool_dict = {}
        subtool_dict['applicationSuite'] = {'name': parent_metadata.name, 'softwareVersion': parent_metadata.softwareVersion, 'identifier': parent_metadata.identifier}
        if not subtool_name in parent_metadata.featureList:
            raise ValueError(f"Cannot create subtool metadata. {subtool_name} is not in featureList of {parent_metadata.name}")
        subtool_dict['name'] = subtool_name
        return cls(**subtool_dict)






class ToolInstanceMetadata(ToolMetadataBase):

    @staticmethod
    def _init_metadata():
        return OrderedDict([
            ('toolName', None),
            ('toolVersion', None),
            ('toolIdentifier', None),
            ('name', None),
            ('version', '0.1.0'),  # Set to something low if not provided.
            ('identifier', None),
            ('description', None),
            ('inputObjects', [IOObjectItem()]),
            ('ouputObjects', [IOObjectItem()]),
        ])
