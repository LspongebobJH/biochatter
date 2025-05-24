import inspect
import json
import re
from types import MappingProxyType, ModuleType
from typing import Any, Callable
import argparse
import importlib

from docstring_parser import parse
from pydantic import Field, create_model, PrivateAttr
from importlib.metadata import version

from langchain_core.utils.pydantic import create_model as create_model_v1
from langchain_core.pydantic_v1 import Field as FieldV1

from datamodel_code_generator import DataModelType, PythonVersion
from datamodel_code_generator.model import get_data_model_types
from datamodel_code_generator.parser.jsonschema import JsonSchemaParser

from .base.agent_abc import BaseAPIModel, BaseAPI

# Jiahang: deprecated in the future, replaced by apis_to_data_models
def generate_pydantic_classes(module: ModuleType) -> list[type[BaseAPIModel]]:
    """Generate Pydantic classes for each callable.

    For each callable (apition/method) in a given module. Extracts parameters
    from docstrings using docstring-parser. Each generated class has fields
    corresponding to the parameters of the apition. If a parameter name
    conflicts with BaseModel attributes, it is aliased.

    Params:
    -------
    module : ModuleType
        The Python module from which to extract apitions and generate models.

    Returns
    -------
    list[Type[BaseModel]]
        A list of Pydantic model classes corresponding to each apition found in
            `module`.

    Notes
    -----
    - For now, all parameter types are set to `Any` to avoid complications with
      complex or external classes that are not easily JSON-serializable.
    - Optional parameters (those with a None default) are represented as
      `Optional[Any]`.
    - Required parameters (no default) use `...` to indicate that the field is
      required.

    """
    base_attributes = set(dir(BaseAPIModel))
    classes_list = []

    for name, api in inspect.getmembers(module, inspect.isapition):
        # Skip private/internal apitions (e.g., _something)
        if name.startswith("_"):
            continue

        # Parse docstring for parameter descriptions
        doc = inspect.getdoc(api) or ""
        parsed_doc = parse(doc)
        doc_params = {p.arg_name: p.description or "No description available." for p in parsed_doc.params}

        sig = inspect.signature(api)
        fields = {}

        for param_name, param in sig.parameters.items():
            # Skip *args and **kwargs for now
            if param_name in ("args", "kwargs"):
                continue

            # Fetch docstring description or fallback
            description = doc_params.get(param_name, "No description available.")

            # Determine default value
            # If no default, we use `...` indicating a required field
            if param.default is not inspect.Parameter.empty:
                default_value = param.default

                # Convert MappingProxyType to a dict for JSON compatibility
                if isinstance(default_value, MappingProxyType):
                    default_value = dict(default_value)

                # Handle non-JSON-compliant float values by converting to string
                if default_value in [float("inf"), float("-inf"), float("nan"), float("-nan")]:
                    default_value = str(default_value)
            else:
                default_value = ...  # No default means required

            # For now, all parameter types are Any
            annotation = Any

            # Append the original annotation as a note in the description if
            # available
            if param.annotation is not inspect.Parameter.empty:
                description += f"\nOriginal type annotation: {param.annotation}"

            # If default_value is None, parameter can be Optional
            # If not required, mark as Optional[Any]
            if default_value is None:
                annotation = Any | None

            # Prepare field kwargs
            field_kwargs = {"description": description, "default": default_value}

            # If field name conflicts with BaseModel attributes, alias it
            field_name = param_name
            if param_name in base_attributes:
                alias_name = param_name + "_param"
                field_kwargs["alias"] = param_name
                field_name = alias_name

            fields[field_name] = (annotation, FieldV1(**field_kwargs))

        # Create the Pydantic model

        tl_parameters_model = create_model_v1(
            name,
            **fields,
            __base__=BaseAPIModel,
        )
        classes_list.append(tl_parameters_model)
    return classes_list

def get_api_path(module: ModuleType, api: Callable) -> str:
    """Get the path of an API.

    This apition takes a module and an API, and returns the path of the API.
    """
    return module.__name__ + '.' + api.__name__

def get_class_name(module: ModuleType, api: Callable) -> str:
    """Get the internal name of an API.

    This apition takes a module and an API, and returns the internal name of the
    API.
    """

    module_name = module.__name__
    api_name = api.__name__

    module_name = ''.join(_name.capitalize() for _name in re.findall(r'[a-zA-Z]+', module_name))
    api_name = ''.join(_name.capitalize() for _name in re.findall(r'[a-zA-Z]+', api_name))
    return module_name + api_name

def get_py_version() -> PythonVersion:
    """Get the Python version.
    """
    from sys import version_info
    py_version = version_info.major, version_info.minor, version_info.micro
    assert py_version[0] == 3, "Python version must be 3.x.x"
    if py_version[1] < 9:
        raise ValueError("Python version must be less than 3.14 and larger than or equal to 3.9.")
    if py_version[1] >= 9 and py_version[1] < 10:
        return PythonVersion.PY_39
    if py_version[1] >= 10 and py_version[1] < 11:
        return PythonVersion.PY_310
    if py_version[1] >= 11 and py_version[1] < 12:
        return PythonVersion.PY_311
    if py_version[1] >= 12 and py_version[1] < 13:
        return PythonVersion.PY_312
    if py_version[1] >= 13 and py_version[1] < 14:
        return PythonVersion.PY_313
    if py_version[1] >= 14:
        raise ValueError("Python version must be less than 3.14 and larger than or equal to 3.9.")

def get_info_import_path(package: ModuleType, object_name: str) -> str:
    package_name = package.__name__
    import_path = f"biochatter.api_agent.python.{package_name}.info_hub.{object_name}"
    return import_path
    
def data_model_to_py(data_model: type[BaseAPIModel], additional_imports: list[str], need_import: bool) -> str:
    """Convert a Pydantic model to a Python code.
    """
    json_schema = json.dumps(data_model.model_json_schema())
    data_model_types = get_data_model_types(
        DataModelType.PydanticV2BaseModel,
        target_python_version=get_py_version()
    )
    parser = JsonSchemaParser(
        json_schema,
        data_model_type=data_model_types.data_model,
        data_model_root_type=data_model_types.root_model,
        data_model_field_type=data_model_types.field_model,
        data_type_manager_type=data_model_types.data_type_manager,
        dump_resolve_reference_action=data_model_types.dump_resolve_reference_action,
        base_class="BaseAPI",
        additional_imports=additional_imports
    )
    codes: str = parser.parse()

    # Jiahang: be noted that following hack could be tricky since it relies on str operations,
    # and thus can be guaranteed on only datamodel_code_generator 0.30.1.

    # remove incorrect import resulted from base_class="BaseAPI" in JsonSchemaParser.
    codes = re.sub(r"^import BaseAPI\s*", "", codes, flags=re.MULTILINE)

    # add private attributes.
    keys = ["_api_name", "_products_original", "_data_name"]
    for key in keys:
        attr = data_model.__private_attributes__[key]
        codes += f"""\n    {key} = PrivateAttr({attr})"""

    # remove model_config as it's already set by base_class.
    codes = re.sub(
        r"model_config\s*=\s*ConfigDict\(\s*.*?\s*\)",
        "",
        codes,
        flags=re.DOTALL
    ).strip()

    # remove import if not needed.
    if not need_import:
        codes = re.sub(r"^\s*(import|from)\s+.*\n?", "", codes, flags=re.MULTILINE)
    
    return codes

def apis_to_data_models(
        api_dict: str, 
        ) -> list[type[BaseAPIModel]]:
    """
    Although we have many string operations like hack in this implementation, all these hacks are bound to
    specific version of datamodel_code_generator and pydantic. They are not bound to any specific package, module
    or API, meaning that they are still generic to any API.
    """
    assert version("datamodel_code_generator") == "0.30.1", \
        "datamodel-code-generator version must be 0.30.1 since some fine-grained operations " \
        "are based on the outputs of this package. Different versions may lead to different outputs " \
        "and thus invalidate those fine-grained operations."
    
    base_attributes = set(dir(BaseAPIModel))
    classes_list = []
    codes_list = []
    api_list = api_dict['api_list']
    package = api_dict['meta']['package']
    module = api_dict['meta']['module']

    _need_import = True

    for _api in api_list:
        api = _api['api']
        assert 'products' in _api and 'data_name' in _api, \
            "configs should contain 'products' and 'data_name'."
        name = api.__name__
        if name.startswith("_"):
            raise Warning(f"apition {name} is private/internal and should not be included in the data model.")

        # Parse docstring for parameter descriptions
        doc = inspect.getdoc(api) or ""
        parsed_doc = parse(doc)
        doc_params = {p.arg_name: p.description or "No description available." for p in parsed_doc.params}

        sig = inspect.signature(api)
        fields = {}

        for param_name, param in sig.parameters.items():
            # Skip *args and **kwargs for now
            if param_name in ("args", "kwargs"):
                continue

            # Fetch docstring description or fallback
            description = doc_params.get(param_name, "No description available.")

            # Determine default value
            # If no default, we use `...` indicating a required field
            if param.default is not inspect.Parameter.empty:
                default_value = param.default

                # Convert MappingProxyType to a dict for JSON compatibility
                if isinstance(default_value, MappingProxyType):
                    default_value = dict(default_value)

                # Handle non-JSON-compliant float values by converting to string
                if default_value in [float("inf"), float("-inf"), float("nan"), float("-nan")]:
                    default_value = str(default_value)
            else:
                default_value = ...  # No default means required

            # For now, all parameter types are Any
            annotation = Any

            # Append the original annotation as a note in the description if
            # available
            if param.annotation is not inspect.Parameter.empty:
                description += f"\nOriginal type annotation: {param.annotation}"

            # If default_value is None, parameter can be Optional
            # If not required, mark as Optional[Any]
            if default_value is None:
                annotation = Any | None

            # Prepare field kwargs
            field_kwargs = {"description": description, "default": default_value}

            # If field name conflicts with BaseModel attributes, alias it
            field_name = param_name
            if param_name in base_attributes:
                alias_name = param_name + "_param"
                field_kwargs["alias"] = param_name
                field_name = alias_name

            fields[field_name] = (annotation, Field(**field_kwargs))

        
        # Create the Pydantic model
        fields['_api_name'] = (str, PrivateAttr(default=get_api_path(module, api)))
        fields['_products_original'] = (str, PrivateAttr(default=_api['products']))
        fields['_data_name'] = (str, PrivateAttr(default=_api['data_name']))

        data_model = create_model(
            get_class_name(module, api),
            __doc__ = doc,
            __base__ = BaseAPI,
            **fields,
        )
        classes_list.append(data_model)

        additional_imports = [
            "biochatter.api_agent.base.agent_abc.BaseAPI",
            "pydantic.PrivateAttr",
        ]
        codes = data_model_to_py(data_model, additional_imports, _need_import)
        codes_list.append(codes)

        # hack. Subsequent codes need no repeated imports. This is important
        # to avoid erros like __future__ import not at the top of the file.
        _need_import = False

    # hack. Add TOOLS_DICT to the end of the file.
    codes = "\n\n".join(codes_list)
    codes += "\n\nTOOLS_DICT = {"
    for data_model in classes_list:
        codes += f"\n    \"{data_model._api_name.default}\": {data_model.__name__},"
    codes += "\n}"

    return classes_list, codes

def get_output_path(package_name: str, api_dict_name: str) -> str:
    return f"biochatter/api_agent/python/{package_name}/{api_dict_name}.py"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--package_name", type=str, required=True)
    parser.add_argument("--api_dict_name", type=str, required=True)
    args = parser.parse_args()

    package_name = args.package_name
    api_dict_name = args.api_dict_name
    output_path = get_output_path(package_name, api_dict_name)

    api_dict = importlib.import_module(f"biochatter.api_agent.python.{package_name}.api_dict")
    api_dict = getattr(api_dict, api_dict_name)

    data_models, codes = apis_to_data_models(api_dict)

    with open(output_path, "w") as f:
        f.write(codes)

    print(f"Data models and codes have been generated and saved to {output_path}.")