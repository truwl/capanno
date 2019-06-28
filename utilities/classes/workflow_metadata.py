#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.

from collections import OrderedDict
from pathlib import Path
from ruamel.yaml import safe_load
from utilities.classes.shared_properties import CodeRepository, WebSite, Person, Publication, Keyword, CallMap

from utilities.classes.metadata_base import MetadataBase

class WorkflowMetadataBase(MetadataBase):
    pass



class WorkflowMetadata(WorkflowMetadataBase):

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
        ('callMap', [CallMap()])
    ])

    @classmethod
    def load_from_file(cls, file_path):
        file_path = Path(file_path)
        with file_path.open('r') as file:
            file_dict = safe_load(file)
        return cls(**file_dict)
