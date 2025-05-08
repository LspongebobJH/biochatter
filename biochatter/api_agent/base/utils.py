from .agent_abc import BaseDependency, BaseAPI, BaseData, BaseKeysInfo, InputAPI, InputDependency   
from copy import deepcopy
import json
import ast

def retrieve_products(src_api: BaseAPI, dep: BaseDependency) -> BaseData:
    """Retrieve products from the src_api._products according to deps.keys_info,
    and assign the retrieved products to self.deps
    
    Note that the deps.keys_info must be a subset of src_api._products.keys_info.
    (Jiahang: this should be asserted in the future)

    For now the retrieval is conducted by removing unneeded objects from src_api._products.data.
    However, I (Jiahang) personally think it's a bit tricky and unsafe.
    We should find a better way to do this.
    """
    def _retrieve_data(src_data: BaseData, src_keys_info: BaseKeysInfo, dst_keys_info: BaseKeysInfo):
        for key, sub_src_keys_info in src_keys_info.keys.items():
            if key not in dst_keys_info.keys:
                if sub_src_keys_info.membership == "item":
                    del src_data[key]
                elif sub_src_keys_info.membership == "attr":
                    delattr(src_data, key)
                elif sub_src_keys_info.membership == "self":
                    # only the first layer can be self
                    # which has been passed
                    raise Exception("Only the first layer can be self membership")
            else:
                # Recursively handle nested structures
                if sub_src_keys_info.membership == "item":
                    sub_src_data = src_data.__getitem__(key)
                    sub_dst_keys_info = dst_keys_info.keys[key]
                elif sub_src_keys_info.membership == "attr":
                    sub_src_data = src_data.__getattribute__(key)
                    sub_dst_keys_info = dst_keys_info.keys[key]
                elif sub_src_keys_info.membership == "self":
                    raise Exception("Only the first layer can be self membership")
                else:
                    raise Exception(f"Invalid membership: {sub_src_keys_info.membership}")
                
                if sub_src_keys_info.keys: # there are objects being included, not just itself.
                    sub_src_data = _retrieve_data(sub_src_data, 
                                                  sub_src_keys_info, sub_dst_keys_info)
               
                if sub_dst_keys_info.membership == "item":
                    src_data.__setitem__(key, sub_src_data)
                elif sub_dst_keys_info.membership == "attr":
                    src_data.__setattr__(key, sub_src_data)
                elif sub_dst_keys_info.membership == "self":
                    raise Exception("Only the first layer can be self membership")
                else:
                    raise Exception(f"Invalid membership: {sub_dst_keys_info.membership}")
        
        return src_data
    
    src_data = deepcopy(src_api._products.data)
    src_keys_info = src_api._products.keys_info
    dst_keys_info = dep.deps.keys_info
    
    # Start recursive copy
    dst_data = _retrieve_data(src_data, 
                              src_keys_info, dst_keys_info)
    
    return dst_data

def get_active_in_apis(deps: list[BaseDependency], api: BaseAPI) -> list[BaseDependency]:
        """Get the active dependencies of the given API."""
        active_deps = []
        for dep in deps:
            if check_active_dep_args(dep.args, api):
                active_deps.append(dep)
        return [dep.u_api for dep in active_deps]

def check_active_dep_args(dep_args: dict, api: BaseAPI) -> bool:
    """Check if the edge required argument is activated by the target api."""
    e_arg_name = list(dep_args.keys())[0] # Jiahang: only one edge required arg being considered for now
    e_arg_val = dep_args[e_arg_name]
    if e_arg_name in api.model_fields.keys() and getattr(api, e_arg_name) == e_arg_val:
        return True
    return False

def _input_product_to_internal_keys_info(
            product: ast.AST, 
            child_name: str | None,
            child_keys_info: BaseKeysInfo | None) -> BaseData:
        """
        The core logic is to convert the products keys info of InputAPI 
        to _products keys info of BaseAPI,
        where the data representations is converted from:
        product = ast.parse("data['d'].e")
        to
        keys_info = BaseKeysInfo(
            membership = "self",
            keys = {
                    "d": BaseKeysInfo(
                        membership = "item",
                        keys = {
                            "e": BaseKeysInfo(
                                membership = "attr",
                            )
                        }
                    )
                }
            )
        
        Jiahang: this function is pretty complicated. Documents need to be carefully revised with
        sufficient examples.
        """
        if child_name is None or child_keys_info is None:
            keys = {}
        else:
            keys = {child_name: child_keys_info}

        if isinstance(product, ast.Name):
            # Base case: just a variable name
            return BaseKeysInfo(
                membership="self",
                keys = keys
            )
        else:
            if isinstance(product, ast.Attribute):
                # Handle attribute access (e.g., .attribute)
                name = product.attr
                keys_info: BaseKeysInfo = BaseKeysInfo(
                    membership="attr",
                    keys = keys
                )
                
            elif isinstance(product, ast.Subscript):
                # Handle subscript access (e.g., ['key'])
                # Jiahang: this assert should be put in field validator of InputAPI
                assert isinstance(product.slice, ast.Constant), "Only constant subscript is supported for now."
                name = product.slice.value
                keys_info: BaseKeysInfo = BaseKeysInfo(
                    membership="item",
                    keys = keys
                )
            else:
                raise Exception(f"Invalid product access: {product}. Only [] and . access are supported.")
            
            parent_product = product.value
            full_keys_info: BaseKeysInfo = \
                _input_product_to_internal_keys_info(
                    parent_product,
                    name,
                    keys_info
                    )
            return full_keys_info
    
def _combine_keys_info(combined: BaseKeysInfo, being_combined: BaseKeysInfo) -> BaseKeysInfo:
    """Combine a list of keys info into a single keys info.
    
    An example of keys_info_list:
    keys_info_list[0] = BaseKeysInfo(
        membership = "self",
        keys = {
                "d": BaseKeysInfo(
                    membership = "item",
                    keys = {
                        "e": BaseKeysInfo(
                            membership = "attr",
                        )
                    }
                )
            }
        )
    
    keys_info_list[1] = BaseKeysInfo(
        membership = "self",
        keys = {
                "d": BaseKeysInfo(
                    membership = "item",
                    keys = {
                        "c": BaseKeysInfo(
                            membership = "attr",
                        )
                    }
                )
            }
        )

    After combining two keys_info, we obtain:
    keys_info = BaseKeysInfo(
        membership = "self",
        keys = {
            "d": BaseKeysInfo(
                membership = "item",
                keys = {
                    "c": BaseKeysInfo(
                        membership = "attr",
                    ),
                    "e": BaseKeysInfo(
                        membership = "attr",
                    )
                }
            }
        )

    Jiahang: this function is pretty complicated. Documents need to be carefully revised with
    sufficient examples.
    """ 
    _combined = BaseKeysInfo(membership=combined.membership)
    _combined.keys = deepcopy(combined.keys)
    
    # Recursively merge any overlapping keys
    for key, value in being_combined.keys.items():
        if key in _combined.keys.keys() and _combined.keys[key].membership == value.membership:
            _combined.keys[key] = _combine_keys_info(_combined.keys[key], value)
        else:
            _combined.keys[key] = value
    return _combined
        
def _str_list_to_keys_info(str_list: list[str]) -> BaseKeysInfo:
    """Convert a string representation of a keys info to a keys info object.
    
    The string representation is a string of the form:
    [
        "data['d'].e",
        "data['d'].f"
    ]

    The keys info is a keys info object of the form:

    BaseKeysInfo(
        membership = "self",
        keys = {
            "d": BaseKeysInfo(
                membership = "item",
                keys = {
                    "e": BaseKeysInfo(
                        membership = "attr",
                    ),
                    "f": BaseKeysInfo(
                        membership = "attr",
                    )
                }
            )
        }
    )

    Jiahang: this function is pretty complicated. Documents need to be carefully revised with
    sufficient examples.
    """
    keys_info_list = []

    for p in str_list:
        p = ast.parse(p).body[0].value
        keys_info: BaseKeysInfo = \
            _input_product_to_internal_keys_info(p, None, None)
        keys_info_list.append(keys_info)

    if len(keys_info_list) == 0:
        result = BaseKeysInfo()
    elif len(keys_info_list) == 1:
        result = keys_info_list[0]
    else:
        result = keys_info_list[0]
        for keys_info in keys_info_list[1:]:
            result = _combine_keys_info(result, keys_info)
    
    return result

def read_apis_from_graph_dict(graph_dict: dict, tools_dict: dict) -> dict[str, BaseAPI]:
    """Jiahang: this function is pretty complicated. Documents need to be carefully revised with
    sufficient examples."""

    apis = {}
    for node in graph_dict["nodes"]:
        input_api = InputAPI.model_validate(node)
        internal_api: BaseAPI = tools_dict[input_api.api]()
        internal_api._api_name = input_api.api

        _products = BaseData()
        _products.keys_info = _str_list_to_keys_info(input_api.products)
        internal_api._products = _products
        apis[internal_api._api_name] = internal_api

    return apis

def read_deps_from_graph_dict(graph_dict: dict) -> dict[str, BaseDependency]:
    """Jiahang: this function is pretty complicated. Documents need to be carefully revised with
    sufficient examples."""

    deps = {}
    for edge in graph_dict["edges"]:
        input_dep = InputDependency.model_validate(edge)
        internal_dep: BaseDependency = BaseDependency(
            u_api_name = input_dep.source,
            v_api_name = input_dep.target,
            args = input_dep.args,
            arg_typs = input_dep.arg_types
        )

        dependencies = BaseData()
        dependencies.keys_info = _str_list_to_keys_info(input_dep.dependencies)
        internal_dep.deps = dependencies
        deps[(internal_dep.u_api_name, internal_dep.v_api_name)] = internal_dep

    return deps

    
    