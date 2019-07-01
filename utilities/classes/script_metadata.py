#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.



from pathlib import Path
from collections import OrderedDict
from ruamel.yaml import safe_load

from utilities.classes.common_functions import is_attr_empty
from utilities.classes.shared_properties import WebSite, CodeRepository, Person, Publication, Keyword, ParentScript, Tool
from utilities.classes.metadata_base import MetadataBase


class ScriptMetadataBase(MetadataBase):
    """Factor stuff out to here if there is more than one ScriptMetadata class."""
    pass


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
        ('parentScripts', ParentScript()),
        ('tools', Tool()),
        ('parentMetadata', None),
        ('_parentMetadata', None),  # Place to store parent ScriptMetadata objects.
    ])

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


