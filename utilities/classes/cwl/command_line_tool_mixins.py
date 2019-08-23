
from cwltool.process import shortname

class CommandLineToolMixin:
    """Mixin methods for working with cwl_classes.CommandLineTool objects. These objects should be preprocessed and
    validated before using these methods"""

    cwl_types = ('null', 'boolean', 'int', 'long', 'float', 'double', 'string', 'File', 'Directory')

    def _get_schema_def_requirement(self):
        """
        Checks for SchemaDefRequirment in requirements, and if it is there returns it.
        :return: SchemaDefRequirement | None
        """
        schema_def_requirement = None
        if self.requirements:
            for requirement in self.requirements:
                try:
                    requirement.__getattribute__('types')
                    schema_def_requirement = requirement
                except AttributeError:
                    pass # Requirement is not SchemaDefRequirement. Ignore it.
        return schema_def_requirement

class CommandInputParameterMixin:

    def _handle_str_input_type(self, _type, schema_def_requirement):
        """
        Should be the terminating function of walking CommandInputParameter.type fields. Return type if type is CWLType,
        otherwise return type from SchemaDefRequirement dict.
        :param _type:
        :param schema_def_dict:
        :return:
        """
        if _type in ('null', 'boolean', 'int', 'long', 'float', 'double', 'string', 'File', 'Directory'):
            _type = _type
        else:
            schema_def_dict = schema_def_requirement._make_schema_def_dict()
            schema_def_name = shortname(_type)
            _type = schema_def_dict[schema_def_name]
        return _type

    def _handle_input_type_field(self, type_field, schema_def_requirement):
        if isinstance(type_field, str):
            input_type = self._handle_str_input_type(type_field, schema_def_requirement)

        elif isinstance(type_field, list):  # got a list of types.
            input_type = []
            for _type in type_field:
                if isinstance(_type, str):
                    input_type.append(self._handle_str_input_type(_type, schema_def_requirement))
                else:  # Need to handle record, enum, and array types here.
                    input_type.append(_type.save())
        else:
            input_type = type_field.save()  # Takes care of CommandInput[Array/Enum/Record]Schmema. Can deal with separately later if we want to.

        return input_type

class SchemaDefRequirementMixin:

    def _make_schema_def_dict(self):
        """
        Make dictionary from SchemaDefRequirement to populate inputs parameters fully.

        :param schema_def_requirement: SchemaDefRequirement object.
        :return: keys are input names, values are input parameter fields to drop into inputs section.
        :rtype: dict
        """
        schema_def_types = self.types  # list of InputRecordSchema | InputArraySchema | InputEnumSchema objects.
        if not isinstance(schema_def_types, list):
            raise NotImplementedError
        schema_def_dict = {}
        for type_def in schema_def_types:
            type_def_name = shortname(type_def.name)
            schema_def_dict[type_def_name] = type_def.save()
            schema_def_dict[type_def_name]['name'] = shortname(schema_def_dict[type_def_name]['name'])
            if type_def.type == 'record':
                pass
            elif type_def.type == 'array':
                pass
            elif type_def.type == 'enum':
                raise NotImplementedError("EnumRecordSchema puts long id's in for symbols. Make sure this is handled.")
            else:
                raise ValueError(f"Unexpected InputSchema type {repr(type_def)}")
        return schema_def_dict