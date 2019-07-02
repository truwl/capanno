# Classes to represent metadata for command line tools.
from pathlib import Path
from collections import OrderedDict
from abc import abstractmethod
import re
from ruamel.yaml import safe_load
from utilities.classes.metadata_base import MetadataBase
from utilities.classes.shared_properties import CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite, \
    IOObjectItem
from utilities.get_metadata_from_biotools import make_tool_metadata_kwargs_from_biotools
from utilities.classes.common_functions import _mk_hashes, NameSoftwareVersionMixin


class ToolMetadataBase(MetadataBase):
    """Factor stuff out to here."""

    @abstractmethod
    def _mk_identifier(self, **kwargs):
        pass

    @abstractmethod
    def _check_identifier(self, identifier):
        pass

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, identifier=None, **kwargs):
        if identifier:
            identifier = self._check_identifier(identifier)
        else:
            identifier = self._mk_identifier(**kwargs)
        self._identifier = identifier


class ToolMetadata(NameSoftwareVersionMixin, ToolMetadataBase):
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


    def _check_identifier(self, identifier):
        if not identifier[:3] == "TL_":
            raise ValueError(f"Tool identifiers must start with 'TL_' you provided {identifier}")
        else:
            hex_pattern = r'[0-9a-f]{6}\.[0-9a-f]{2}$'
            match_obj = re.match(hex_pattern, identifier[3:])
            if not match_obj:
                raise ValueError(f"Tool identifier not formatted correctly: {identifier}")

        return identifier

    def _mk_identifier(self, start=0):
        if not (self.name and self.softwareVersion):
            raise ValueError(f"Name and softwareVersion must be provided to make an identifier.")
        name_hash, version_hash = _mk_hashes(self.name, self.softwareVersion)
        identifier = f"TL_{name_hash[start:start + 6]}.{version_hash[:2]}"
        return identifier

    # Class factory methods

    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict)

    @classmethod
    def create_from_biotools(cls, biotools_id, softwareVersion, version='0.1.1'):
        kwargs = make_tool_metadata_kwargs_from_biotools(biotools_id)
        kwargs['softwareVersion'] = softwareVersion
        kwargs['version'] = version
        return cls(**kwargs)


    def mk_instance(self):
        tool_inst = ToolInstanceMetadata(toolName=self.name, toolVersion=self.version, toolIdentifier=self.identifier)
        return tool_inst

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


    def _check_identifier(self, identifier):
        if not identifier[:3] == "TL_":
            raise ValueError(f"Tool identifiers must start with 'TL_' you provided {identifier}")
        else:
            hex_pattern = r'[0-9a-f]{6}\.[0-9a-f]{2}$'
            match_obj = re.match(hex_pattern, identifier[3:])
            if not match_obj:
                raise ValueError(f"Tool identifier not formatted correctly: {identifier}")

        return identifier

    def _mk_identifier(self, start=0):
        if not (self.name and self.softwareVersion):
            raise ValueError(f"Name and softwareVersion must be provided to make an identifier.")
        name_hash, version_hash = _mk_hashes(self.name, self.softwareVersion)
        identifier = f"TL_{name_hash[start:start + 6]}.{version_hash[:2]}"
        return identifier

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
    def create_from_biotools(cls, biotools_id, softwareVersion, subtools, version='0.1.1'):
        kwargs = make_tool_metadata_kwargs_from_biotools(biotools_id)
        kwargs['featureList'] = list(subtools)
        kwargs['softwareVersion'] = softwareVersion
        kwargs['version'] = version
        return cls(**kwargs)


class SubtoolMetadata(MetadataBase):

    @staticmethod
    def _init_metadata():
        return OrderedDict([
        ('applicationSuite', ApplicationSuite()),
        ('name', None),
        ('version', '0.1'),
        ('description', None),
        ('keywords', Keyword()),
        ('alternateName', None),
        ('_parentMetadata', None),
    ])

    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict)


    def _mk_identifier(self):
        raise NotImplementedError

    def _check_identifier(self, identifier):
        raise NotImplementedError

    @classmethod
    def initialize_from_parent(cls, parent_metadata, subtool_name):
        subtool_dict = {}
        subtool_dict['applicationSuite'] = {'name': parent_metadata.name, 'softwareVersion': parent_metadata.softwareVersion, 'identifier': parent_metadata.identifier}
        if not subtool_name in parent_metadata.featureList:
            raise ValueError(f"Cannot create subtool metadata. {subtool_name} is not in featureList of {parent_metadata.name}")
        subtool_dict['name'] = subtool_name
        return cls(**subtool_dict)

    def mk_completed_file(self):
        # Will contain attributes defined in subtool AND from parent.
        raise NotImplementedError


    def mk_instance(self):
        raise NotImplementedError


class ToolInstanceMetadata(ToolMetadataBase):

    # instances of ToolInstanceMetadata should not be instantiated directly, but instantiated from ToolMetadata and SubtoolMetadata.
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
            ('_tool_metadata', None)
        ])

