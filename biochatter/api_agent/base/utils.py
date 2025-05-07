from .agent_abc import BaseDependency, BaseAPI, BaseData, BaseKeysInfo
from copy import deepcopy

def retrieve_products(src_api: BaseAPI, dep: BaseDependency) -> BaseData:
    """Retrieve products from the src_api._products according to deps.keys_info,
    and assign the retrieved products to self.deps
    
    Note that the deps.keys_info must be a subset of src_api._products.keys_info

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
                
                if sub_src_keys_info.keys is not None: # there are objects being included, not just itself.
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
