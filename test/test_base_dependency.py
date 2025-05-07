from unittest.mock import MagicMock

from biochatter.api_agent.base.agent_abc import BaseAPI, BaseDependency, BaseData, BaseKeysInfo
from biochatter.api_agent.base.utils import retrieve_products
def test_data():
    
    class TestData:
        def __init__(self, initialize = True):
            if initialize:
                self.obj1 = MagicMock()
                self.obj2 = {
                    "key1": MagicMock(),
                    "key2": MagicMock(),
                }

                test_obj = MagicMock()
                test_obj.key4 = MagicMock()
                self.data_dict = {
                    "key3": test_obj,
                }

        def __getitem__(self, name):
            return self.data_dict[name]
        
        def __setitem__(self, name, value):
            self.data_dict[name] = value

    keys_info = {
        "membership": "self",
        "keys": {
            "obj1": {
                "membership": "attr",
            },
            "obj2": {
                "membership": "attr",
                "keys": {
                    "key1": {
                        "membership": "item",
                    },
                    "key2": {
                        "membership": "item",
                    },
                },
            },
            "key3": {
                "membership": "item",
                "keys": {
                    "key4": {
                        "membership": "attr",
                    },
                },
            },
        }
    }

    keys_info = BaseKeysInfo.model_validate(keys_info)

    data = BaseData(
        data=TestData(),
        keys_info=keys_info
    )

    src_api = BaseAPI()
    src_api._products = data

    assert data.data.obj1
    assert data.data.obj2["key1"]
    assert data.data.obj2["key2"]
    assert data.data["key3"].key4

    dep_keys_info = {
        "membership": "self",
        "keys": {
            "obj1": {
                "membership": "attr"
            },
            "obj2": {
                "membership": "attr",
                "keys": {
                    "key1": {
                        "membership": "item"
                    }
                }
            },
            "key3": {
                "membership": "item",
                "keys": {
                    "key4": {
                        "membership": "attr"
                    }
                }
            }
        }
    }

    dep_keys_info = BaseKeysInfo.model_validate(dep_keys_info, strict=True)

    dep_data = BaseData(
        data=TestData(initialize=False),
        keys_info=dep_keys_info
    )


    dependency = BaseDependency(
        u_api_name="src",
        v_api_name="dst",
        args={},
        arg_typs={},
        deps=dep_data,
    )

    result = retrieve_products(src_api, dependency)

    assert result.obj1 == src_api._products.data.obj1
    assert result.obj2["key1"] == src_api._products.data.obj2["key1"]
    assert "key2" not in result.obj2.keys()
    assert result["key3"].key4 == src_api._products.data["key3"].key4


if __name__ == "__main__":
    test_data()

    pass