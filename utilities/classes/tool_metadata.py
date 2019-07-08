# Classes to represent metadata for command line tools.
from pathlib import Path
from collections import OrderedDict
from abc import abstractmethod
import re
from ruamel.yaml import safe_load
from utilities.classes.metadata_base import MetadataBase
from utilities.classes.shared_properties import CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite, \
    IOObjectItem
from utilities.helpers.get_metadata_from_biotools import make_tool_metadata_kwargs_from_biotools
from utilities.classes.common_functions import _mk_hashes, NameSoftwareVersionMixin, is_attr_empty


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
        raise NotImplementedError


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

    def make_subtool_metadata(self, subtool_name, parent_metadata_path):
        if not self.featureList:
            raise ValueError(f"Cannot create subtool. featureList of {self.name} is not populated.")
        if subtool_name not in self.featureList:
            raise ValueError(f"{subtool_name} must be in the parent featureList")
        subtool_metadata = SubtoolMetadata(name=subtool_name, _parentMetadata=self, parentMetadata=parent_metadata_path)
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
            ('identifier', None),
            ('description', None),
            ('keywords', Keyword()),
            ('alternateName', None),
            ('parentMetadata', None),  # relative path to parentMetadata
            ('_parentMetadata', None),  # ParentMetadata instance. Can be loaded from parentMetadata or set directly.
            ('_primary_file_attrs', None),
        ])

    def __init__(self, file_path=None, **kwargs):
        # file_path only used if class is initiated using 'load_from_file'. file_path is path that SubtoolMetadata is loaded from.
        self._parentMetadata = kwargs.get('_parentMetadata')
        if not self._parentMetadata:
            self.parentMetadata = kwargs['parentMetadata']  # must have a path if it isn't set directly.
            self._load_parent_metadata(file_path)  # sets self._parentMetadata
        self._primary_file_attrs = []
        for k, value in kwargs.items():
            if value:
                self._primary_file_attrs.append(k)  # keep track of kwargs supplied.
        self._load_attrs_from_parent()
        super().__init__(**kwargs)


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not value in self._parentMetadata.featureList:
            raise ValueError(f"{value} is not in {self._parentMetadata.name} metadata featureList")
        self._name = value


    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict, file_path=file_path)

    def _load_parent_metadata(self, subtool_metadata_file_path):
        # subtool_metadata_file_path of subtool metadata. ParentMetadata path should be relative to this.
        dir_name = subtool_metadata_file_path.parent
        full_path = dir_name / self.parentMetadata
        with full_path.resolve().open('r') as f:
            parent_metadata_dict = safe_load(f)
        self._parentMetadata = ParentToolMetadata(**parent_metadata_dict)

    def _load_attrs_from_parent(self):
        # initialize everything from parent. Will be overwritten anything supplied in kwargs.
        parent_meta = self._parentMetadata
        self.applicationSuite = ApplicationSuite(name=parent_meta.name, softwareVersion=parent_meta.softwareVersion,
                                                 identifier=parent_meta.identifier)
        # self.identifier = self._mk_identifier()
        self.keywords = parent_meta.keywords
        return

    def _mk_identifier(self):
        identifier_str, version_str = self._parentMetadata.identifier.split('.', 1)
        subtool_hash = _mk_hashes(self.name)
        identifier = f"{identifier_str}_{subtool_hash}.{version_str}"
        return identifier

    def _check_identifier(self, identifier):
        parent_identifier = self._parentMetadata.identifier
        if not identifier.startswith(parent_identifier[:10]):
            raise ValueError(f"Subtool identifier {identifier} does not properly correspond to parent identifier {parent_identifier}")
        if not identifier.endswith(parent_identifier[-3:]):  # should be '.xx'
            raise ValueError(
                f"Subtool identifier {identifier} does not properly correspond to parent identifier {parent_identifier}")
        hex_pattern = r'[0-9a-f]{6}_[0-9a-f]{2}\.[0-9a-f]{2}$'
        match_obj = re.match(hex_pattern, identifier[3:])
        if not match_obj:
            raise ValueError(f"Tool identifier not formatted correctly: {identifier}")

        return identifier

    @classmethod
    def initialize_from_parent(cls, parent_metadata, subtool_name):
        subtool_dict = {}
        subtool_dict['applicationSuite'] = {'name': parent_metadata.name,
                                            'softwareVersion': parent_metadata.softwareVersion,
                                            'identifier': parent_metadata.identifier}
        if not subtool_name in parent_metadata.featureList:
            raise ValueError(
                f"Cannot create subtool metadata. {subtool_name} is not in featureList of {parent_metadata.name}")
        subtool_dict['name'] = subtool_name
        return cls(**subtool_dict)

    def mk_completed_file(self):
        raise NotImplementedError

    def mk_instance(self):
        raise NotImplementedError

    def mk_file(self, file_path, keys=None):
        parent_path = self.parentMetadata
        subtool_path = file_path
        rel_path = parent_path.relative_to(subtool_path.parents[1])
        rel_path_str = '../' + str(rel_path)
        self.parentMetadata = rel_path_str
        super().mk_file(file_path)
