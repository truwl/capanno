
from abc import abstractmethod
from pathlib import Path
import re
from collections import OrderedDict
from ruamel.yaml import safe_load
from utilities.classes.common_functions import is_attr_empty, NameSoftwareVersionMixin, _mk_hashes
from utilities.classes.shared_properties import WebSite, CodeRepository, Person, Publication, Keyword, ParentScript, Tool
from utilities.classes.metadata_base import MetadataBase


class ScriptMetadataBase(MetadataBase):
    """Factor stuff out to here if there is more than one ScriptMetadata class."""
    pass




class ScriptMetadata(NameSoftwareVersionMixin, ScriptMetadataBase):

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
        ('parentScripts', ParentScript()),
        ('tools', Tool()),
        ('parentMetadata', None),
        ('_parentMetadata', None),  # Place to store parent ScriptMetadata objects.
    ])

    def _check_identifier(self, identifier):
        if not identifier[:3] == "ST_":
            raise ValueError(f"Script identifiers must start with 'ST_' you provided {identifier}")
        else:
            hex_pattern = r'[0-9a-f]{6}\.[0-9a-f]{2}$'
            match_obj = re.match(hex_pattern, identifier[3:])
            if not match_obj:
                raise ValueError(f"Script identifier not formatted correctly: {identifier}")

        return identifier

    def _mk_identifier(self, start=0):
        if not (self.name and self.softwareVersion):
            raise ValueError(f"Name and softwareVersion must be provided to make an identifier.")
        name_hash, version_hash = _mk_hashes(self.name, self.softwareVersion)
        identifier = f"ST_{name_hash[start:start + 6]}.{version_hash[:2]}"
        return identifier

    def _load_common_metadata(self, file_path):
        # Start with returning a list of dicts.
        if self.parentMetadata:
            common_meta_list = []  # List of parent ScriptMetadata objects.
            dirname = file_path.parents[0]
            for rel_path in self.parentMetadata:
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
        new_instance._parentMetadata = new_instance._load_common_metadata(file_path)
        return new_instance


    def _update_attributes(self, update_instance):
        """
        Update self with attribute values in update_instance for attributes in self that have not been set.
        :param update_instance:
        :return:
        """

        for attribute_name in self._init_metadata().keys():
            if is_attr_empty(getattr(self, attribute_name)):
                update_value = getattr(update_instance, attribute_name)
                if is_attr_empty(update_value):
                    continue
                else:
                    setattr(self, attribute_name, update_value)
            else:
                continue
        return


    def mk_completed_file(self, file_path):
        # substitute in parent metadata fields for fields not specified by the script's own metadata.
        for parent_metadata in self._parentMetadata:
            self._update_attributes(parent_metadata)
            MetadataBase.mk_file(self, file_path)
        return

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


class CommonScriptMetadata(ScriptMetadataBase):
    @staticmethod
    def _init_metadata():
        # Only attributes which can be common to multiple scripts.
        return OrderedDict([
            ('softwareVersion', None),
            ('description', None),
            ('WebSite', [WebSite()]),
            ('codeRepository', CodeRepository()),
            ('license', None),
            ('contactPoint', [Person()]),
            ('publication', [Publication()]),
            ('keywords', [Keyword()]),
            ('creator', [Person()]),
            ('programmingLanguage', None),
            ('datePublished', None),
            ('parentScripts', ParentScript()),
            ('tools', Tool()),
            ])
