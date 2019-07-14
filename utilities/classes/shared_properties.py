#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.md', which is part of this source code package.

from abc import abstractmethod, ABC
from urllib.parse import urlparse
from ruamel.yaml.comments import CommentedMap

class AttributeBase(ABC):

    def __init__(self, *args, **kwargs):
        __slots__ = list(self.attrs)

    def is_empty(self):
        _is_empty = True
        for attribute in self.attrs:
            if isinstance(attribute, list):
                for item in attribute:
                    if getattr(item, 'attrs'):
                        if not item.is_empty():
                            _is_empty = False
                            break
                    elif item:
                        _is_empty = False
                        break
            elif getattr(self, attribute):
                _is_empty = False
                break
            else:
                pass
        return _is_empty

    @staticmethod
    @abstractmethod
    def _attrs():
        return frozenset([])


    @property
    def attrs(self):
        return self._attrs()

    def dump(self):
        map_object = CommentedMap()
        for attribute in self.attrs:
            map_object[attribute] = getattr(self, attribute)
        return map_object


class CodeRepository(AttributeBase):
    def __init__(self, name=None, URL=None):
        super().__init__()
        self._name = name
        self._URL = URL
        return

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name
        return

    @property
    def URL(self):
        return self._URL

    @URL.setter
    def URL(self, new_URL):
        if new_URL:
            valid_schemes = ['https', 'http', 'git']
            parse_result = urlparse(new_URL)
            if parse_result.scheme not in valid_schemes:
                raise ValueError(f"URL scheme should be in {valid_schemes}")
        else:
            new_URL = None
        self._URL = new_URL
        return


    @staticmethod
    def _attrs():
        return frozenset(['name', 'URL'])


class WebSite(AttributeBase):
    def __init__(self, name=None, description=None, URL=None):
        super().__init__()
        self._name = name
        self._description = description
        self._URL = URL

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, new_description):
        self._description = new_description

    @property
    def URL(self):
        return self._URL

    @URL.setter
    def URL(self, new_URL):
        if new_URL:
            valid_schemes = ['https', 'http']
            parse_result = urlparse(new_URL)
            if parse_result.scheme not in valid_schemes:
                raise ValueError(f"URL scheme should be in {valid_schemes}")
        else:
            new_URL = None
        self._URL = new_URL
        return

    @staticmethod
    def _attrs():
        return frozenset(['name', 'description', 'URL'])


class Publication(AttributeBase):
    def __init__(self, identifier=None, headline=None):
        super().__init__()
        self._identifier = identifier
        self._headline = headline
        return

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier
        return

    @property
    def headline(self):
        return self._headline

    @headline.setter
    def headline(self, new_headline):
        self._headline = new_headline
        return

    @staticmethod
    def _attrs():
        return frozenset(['identifier', 'headline'])


class Person(AttributeBase):
    def __init__(self, name=None, email=None, identifier=None):
        super().__init__()
        self._name = name
        self._email = email
        self._identifier = identifier

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def email(self):
        return self._email

    @email.setter
    def email(self, new_email):
        self._email = new_email
        return

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier

    @staticmethod
    def _attrs():
        return frozenset(['name', 'email', 'identifier'])


class Keyword(AttributeBase):
    def __init__(self, *args, **kwargs):  # Need to initialize off of *args and **kwargs to handle both forms.
        super().__init__()
        args_len = len(args)
        if args_len == 0:
            self._uri = kwargs.get('uri', None)
        elif args_len == 1:
            self._uri = args[0]
        else:
            raise ValueError(f"Expected only one argument for uri for keyword. Got {args}")
        self._name = kwargs.get('name', None)
        self._category = kwargs.get('category', None)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def category(self):
        return self._category

    @category.setter
    def category(self, new_category):
        category_values = (None, 'topic', 'operation')
        if new_category not in category_values:
            raise ValueError(f"{new_category} is not a valid category for a keyword. Must be one of {category_values}")
        self._category = new_category

    @property
    def uri(self):
        return self._uri

    @uri.setter
    def uri(self, new_uri):
        self._uri = new_uri

    def dump(self):
        if self.uri:
            return self.uri
        else:
            keyword = CommentedMap([('name', self.name), ('category', self.category)])
            return keyword

    @staticmethod
    def _attrs():
        return frozenset(['name', 'category', 'uri'])

class ApplicationSuite(AttributeBase):

    def __init__(self, name=None, softwareVersion=None, identifier=None):
        super().__init__()
        self._name = name
        self._softwareVersion = softwareVersion
        self._identifier = identifier

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def softwareVersion(self):
        return self._softwareVersion

    @softwareVersion.setter
    def softwareVersion(self, new_softwareVersion):
        self._softwareVersion = new_softwareVersion

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier


    @staticmethod
    def _attrs():
        return frozenset(['name', 'softwareVersion', 'identifier'])


class ParentScript(AttributeBase):

    def __init__(self, name=None, softwareVersion=None, identifier=None):
        super().__init__()
        self._name = name
        self._softwareVersion = softwareVersion
        self._identifier = identifier

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name


    @property
    def softwareVersion(self):
        return self._softwareVersion

    @softwareVersion.setter
    def softwareVersion(self, new_softwareVersion):
        self._softwareVersion = new_softwareVersion

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier

    @staticmethod
    def _attrs():
        return frozenset(['name','softwareVersion', 'identifier'])


class Tool(AttributeBase):

    def __init__(self, name=None, softwareVersion=None, identifier=None, alternateName=None):
        super().__init__()
        self._name = name
        self._alternateName = alternateName
        self._softwareVersion = softwareVersion
        self._identifier = identifier

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def alternateName(self):
        return self._alternateName

    @alternateName.setter
    def alternateName(self, value):
        self._alternateName = value

    @property
    def softwareVersion(self):
        return self._softwareVersion

    @softwareVersion.setter
    def softwareVersion(self, new_softwareVersion):
        self._softwareVersion = new_softwareVersion

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier


    @staticmethod
    def _attrs():
        return frozenset(['name', 'alternateName', 'softwareVersion', 'identifier'])


class CallMap(AttributeBase):

    def __init__(self, id_=None, identifier=None):
        super().__init__()

        self._id = id_
        self._identifier = identifier

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        self._id = new_id

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier

    @staticmethod
    def _attrs():
        return frozenset(['id', 'identifier'])


class IOObject(AttributeBase):

    def __init__(self, identifier=None, path=None):
        super().__init__()
        self._identifier = identifier
        self._path = path

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._identifier = new_identifier

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, new_path):
        self._path = new_path


    @staticmethod
    def _attrs():
        return frozenset(['identifier', 'path'])

class IOObjectItem(AttributeBase):
    def __init__(self, id_=None, io_object=IOObject()):
        super().__init__()

        self._id = id_
        self._io_object = io_object

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        self._id = new_id

    @property
    def identifier(self):
        return self._io_object._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        self._io_object._identifier = new_identifier

    @property
    def path(self):
        return self._io_object._path

    @path.setter
    def path(self, new_path):
        self._io_object._path = new_path

    @property
    def io_object(self):
        return self._io_object

    @staticmethod
    def _attrs():
        return frozenset(['id', 'io_object'])

    def dump(self):
        map_object = CommentedMap()
        map_object['id'] = getattr(self, 'id')
        map_object.update(self.io_object.dump())
        return map_object

class IOArrayItem(AttributeBase):
    def __init__(self, id_, objects=None):
        super().__init__()
        self._id = id_
        self._objects = objects if objects else [IOObject()]

    @staticmethod
    def _attrs():
        return frozenset(['id', 'objects'])




