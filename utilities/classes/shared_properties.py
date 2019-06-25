#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.md', which is part of this source code package.

from urllib.parse import urlparse
from ruamel.yaml.comments import CommentedMap


class CodeRepository:
    def __init__(self, name=None, URL=None):
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

    def dump(self):
        code_repo = CommentedMap([('name', self.name), ('URL', self.URL)])
        return code_repo


class WebSite:
    def __init__(self, name=None, description=None, URL=None):
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

    def dump(self):
        web_site = CommentedMap([('name', self.name), ('description', self.description), ('URL', self.URL)])
        return web_site


class Publication:
    def __init__(self, identifier=None, headline=None):
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

    def dump(self):
        publication = CommentedMap([('identifier', self.identifier), ('headline', self.headline)])
        return publication


class Person:
    def __init__(self, name=None, email=None, identifier=None):
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

    def dump(self):
        person = CommentedMap([('name', self.name), ('email', self.email), ('identifier', self.identifier)])
        return person


class Keyword:
    def __init__(self, *args, **kwargs):  # Might need to initialize off of *args and **kwargs to handle both forms.
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


class ApplicationSuite:

    def __init__(self, name=None, softwareVersion=None, identifier=None):
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

    def dump(self):
        application_suite = CommentedMap(
            [('name', self.name), ('softwareVersion', self.softwareVersion), ('identifier', self.identifier)])
        return application_suite
