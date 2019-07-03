#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.


from abc import ABC, abstractmethod
import logging
from pathlib import Path
import semantic_version
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml import YAML
from utilities.classes.shared_properties import CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite, ParentScript, Tool, IOObjectItem

object_attributes = (CodeRepository, Person, Publication, WebSite, Keyword, ApplicationSuite, ParentScript, Tool, IOObjectItem)

class MetadataBase(ABC):
    """Factor stuff out to here."""

    @staticmethod
    @abstractmethod
    def _init_metadata():
        return dict([('name', None), ('version', None)])

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
        init_metadata = self._init_metadata()

        for k, v in kwargs.items():
            if not k in init_metadata:
                raise AttributeError(f"{k} is not a valid key for {type(self)}")

        for k, v in init_metadata.items():
            if k in kwargs:
                setattr(self, k, kwargs[k])  # Highest priority.
            else:
                try:
                    if getattr(self, k):  # value has already been set by derived class __init__. Second highest priority.
                        continue
                    else:
                        raise NotImplementedError(f"Figure out what's happening here and fix it.")
                except AttributeError:
                    setattr(self, k, v)  # Set to default value provided in self._init_metadata. Last resort.
        return

    def mk_file(self, file_path, keys=None):
        file_path = Path(file_path)
        meta_map = CommentedMap()
        if not keys:
            keys = self._get_metafile_keys()
        for key in keys:
            if key.startswith('_'):
                continue
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