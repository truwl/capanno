#
# * This file is subject to the terms and conditions defined in
# * file 'LICENSE.txt', which is part of this source code package.

from utilities.classes.metadata_base import object_attributes
from hashlib import md5
import uuid
import base64


def _mk_hashes(arg1, *args):
    arg1 = str(arg1)
    hashes = [md5(arg1.encode(encoding='utf-8')).hexdigest()]
    for arg in args:
        arg = str(arg)
        hashes.append(md5(arg.encode(encoding='utf-8')).hexdigest())
    return hashes

def mk_tool_identifier(name, version, start=0):
    name_hash, version_hash = _mk_hashes(name, version)
    identifier = f"TL_{name_hash[start:start+6]}.{version_hash[:2]}"
    return identifier

def mk_tool_instance_identifier(tool_identifier):
    uuid_string = uuid.uuid4().hex[:4]
    return f"{tool_identifier}.{uuid_string}"

def mk_subtool_identifier():
    raise NotImplementedError

def mk_script_identifier():
    raise NotImplementedError

def mk_workflow_identifier():
    raise NotImplementedError


def is_attr_empty(attribute):
    # Need to check for lists.
    if isinstance(attribute, list):
        is_empty = True
        for item in attribute:
            if not is_attr_empty(item):
                is_empty = False
                break
    elif isinstance(attribute, object_attributes):
        is_empty = attribute.is_empty()
    else:
        if attribute:
            is_empty = False
        else:
            is_empty = True
    return is_empty