import unittest
from unittest.mock import MagicMock
from copy import deepcopy

from biochatter.api_agent.base.agent_abc import BaseAPI, BaseDependency, BaseData, BaseKeysInfo
from biochatter.api_agent.base.utils import (
    retrieve_products, 
    read_apis_from_graph_dict, 
    read_deps_from_graph_dict, 
    aggregate_deps
)
from biochatter.api_agent.python.scanpy import SCANPY_TOOLS_DICT
from biochatter.api_agent.dep_graph.graph import DependencyGraph

class BaseTestData:
    def __init__(self):
        pass

    def __getitem__(self, name):
        return self.data_dict[name]
    
    def __setitem__(self, name, value):
        self.data_dict[name] = value
        
class TestData(BaseTestData):
    def __init__(self):
        super().__init__()
        self.obj1 = "obj1"
        self.obj2 = {
            "key1": "obj2_key1",
            "key2": "obj2_key2",
        }

        self.key_3_obj = MagicMock()
        self.key_3_obj.obj4 = "key3_obj4"
        self.data_dict = {
            "key3": self.key_3_obj,
        }

        self.keys_info = BaseKeysInfo.model_validate(
            {
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
                            "obj4": {
                                "membership": "attr",
                            },
                        },
                    },
                }
            }
        )

class TestDataSrc1(BaseTestData):
    def __init__(self):
        super().__init__()
        self.obj1 = "obj1"
        self.obj2 = {
            "key1": "obj2_key1",
            "key2": "obj2_key2",
            "key3": "obj2_key3",
            "key4": "obj2_key4",
        }
        self.obj4 = {
            "key2": "obj4_key2",
        }
        self.key_3_obj = MagicMock()
        self.key_3_obj.obj6 = "key3_obj6"
        self.data_dict = {
            "key3": self.key_3_obj,
        }

        self.keys_info = BaseKeysInfo.model_validate(
            {
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
                            "key3": {
                                "membership": "item",
                            },
                            "key4": {
                                "membership": "item",
                            },
                        },
                    },
                    "obj4": {
                        "membership": "attr",
                        "keys": {
                            "key2": {
                                "membership": "item",
                            },
                        },
                    },
                    "key3": {
                        "membership": "item",
                        "keys": {
                            "obj6": {
                                "membership": "attr",
                            },
                        },
                    },
                },
            },
        )

class TestDataSrc2(BaseTestData):
    def __init__(self):
        super().__init__()
        self.obj1 = "obj1"
        self.obj2 = {
            "key1": "obj2_key1",
        }
        self.obj3 = "obj3"
        self.obj4 = {
            "key1": "obj4_key1",
            "key2": "obj4_key2",
        }
        self.key_3_obj = MagicMock()
        self.key_3_obj.obj4 = "key3_obj4"
        self.key_3_obj.obj5 = "key3_obj5"
        self.data_dict = {
            "key3": self.key_3_obj,
        }

        self.keys_info = BaseKeysInfo.model_validate(
            {
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
                        },
                    },
                    "obj3": {
                        "membership": "attr",
                    },
                    "obj4": {
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
                            "obj4": {
                                "membership": "attr",
                            },
                            "obj5": {
                                "membership": "attr",
                            },
                        },
                    },
                },
            }
        )

class TestDataGT(BaseTestData):
    def __init__(self):
        super().__init__()
        self.obj1 = "obj1"
        self.obj2 = {
            "key1": "obj2_key1",
            "key2": "obj2_key2",
            "key3": "obj2_key3",
            "key4": "obj2_key4",
        }
        self.obj3 = "obj3"
        self.obj4 = {
            "key1": "obj4_key1",
            "key2": "obj4_key2",
        }

        self.key_3_obj = MagicMock()
        self.key_3_obj.obj4 = "key3_obj4"
        self.key_3_obj.obj5 = "key3_obj5"
        self.key_3_obj.obj6 = "key3_obj6"
        self.data_dict = {
            "key3": self.key_3_obj,
        }

        self.keys_info = BaseKeysInfo.model_validate(
                {
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
                                "key3": {
                                    "membership": "item",
                                },
                                "key4": {
                                    "membership": "item",
                                },
                            },
                        },
                        "obj3": {
                            "membership": "attr",
                        },
                        "obj4": {
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
                                "obj4": {
                                    "membership": "attr",
                                },
                                "obj5": {
                                    "membership": "attr",
                                },
                                "obj6": {
                                    "membership": "attr",
                                },
                            },
                        },
                    },
                },
            )

    def __eq__(self, other):
        assert self.obj1 == other.obj1
        assert self.obj2 == other.obj2
        assert self.obj3 == other.obj3
        assert self.obj4 == other.obj4
        assert self.data_dict['key3'].obj4 == other.data_dict['key3'].obj4
        assert self.data_dict['key3'].obj5 == other.data_dict['key3'].obj5
        assert self.data_dict['key3'].obj6 == other.data_dict['key3'].obj6
        return True


class TestDepGraph(unittest.TestCase):
    def test_construct_dep_graph(self):
        mock_graph_dict = {
            "directed": True,
            "multigraph": False,
            "graph": {},
            "nodes": [
                {
                    "api": "sc.pp.neighbors",
                    "products": [
                        "data.obsm['X_pca']",
                        "data.varm['PCs']",
                        "data.uns['pca']",
                        "data.obsp['distances']",
                        "data.obsp['connectivities']"
                    ],
                    "id": "sc.pp.neighbors"
                },
                {
                    "api": "sc.tl.umap",
                    "products": [
                        "data.obsm['X_umap']",
                        "data.uns['umap']"
                    ],
                    "id": "sc.tl.umap"
                },
                {
                    "api": "sc.tl.leiden",
                    "products": [
                        "data.obs['leiden']",
                        "data.uns['leiden']"
                    ],
                    "id": "sc.tl.leiden"
                },
                {
                    "api": "sc.pl.umap",
                    "products": [],
                    "id": "sc.pl.umap"
                },
                {   
                    "api": "root",
                    "products": [
                        "data.X"
                    ],
                    "id": "root"
                }
            ],
            "edges": [
                {
                    "dependencies": [
                        "data.X"
                    ],
                    "source": "root",
                    "target": "sc.pp.neighbors",
                    "args": {
                        "data": "data"
                    },
                    "arg_types": {
                        "data": "object"
                    }
                },
                {
                    "dependencies": [
                        "data.obsp['distances']"
                    ],
                    "source": "sc.pp.neighbors",
                    "target": "sc.tl.leiden",
                    "args": {
                        "data": "data"
                    },
                    "arg_types": {
                        "data": "object"
                    }
                },
                {
                    "dependencies": [
                        "data.uns['neighbors']",
                        "data.obsp['connectivities']"
                    ],
                    "source": "sc.pp.neighbors",
                    "target": "sc.tl.umap",
                    "args": {
                        "data": "data"
                    },
                    "arg_types": {
                        "data": "object"
                    }
                },
                {
                    "dependencies": [
                        "data.obsm['X_umap']"
                    ],
                    "source": "sc.tl.umap",
                    "target": "sc.pl.umap",
                    "args": {
                        "data": "data"
                    },
                    "arg_types": {
                        "data": "object"
                    }
                },
                {
                    "dependencies": [
                        "data.obs['leiden']"
                    ],
                    "source": "sc.tl.leiden",
                    "target": "sc.pl.umap",
                    "args": {
                        "color": "leiden"
                    },
                    "arg_types": {
                        "color": "str"
                    }
                }
            ]
        }

        graph = DependencyGraph(mock_graph_dict, SCANPY_TOOLS_DICT)

        gt_apis_dict = {
            "sc.pp.neighbors": SCANPY_TOOLS_DICT["sc.pp.neighbors"](),
            "sc.tl.umap": SCANPY_TOOLS_DICT["sc.tl.umap"](),
            "sc.tl.leiden": SCANPY_TOOLS_DICT["sc.tl.leiden"](),
            "sc.pl.umap": SCANPY_TOOLS_DICT["sc.pl.umap"](),
            "root": SCANPY_TOOLS_DICT["root"](),
        }

        gt_apis_dict["sc.pp.neighbors"]._api_name = "sc.pp.neighbors"
        gt_apis_dict["sc.tl.umap"]._api_name = "sc.tl.umap"
        gt_apis_dict["sc.tl.leiden"]._api_name = "sc.tl.leiden"
        gt_apis_dict["sc.pl.umap"]._api_name = "sc.pl.umap"
        gt_apis_dict["root"]._api_name = "root"

        gt_apis_dict["sc.pp.neighbors"]._products = BaseData()
        gt_apis_dict["sc.pp.neighbors"]._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obsm": {
                        "membership": "attr",
                        "keys": {
                            "X_pca": {
                                "membership": "item",
                            }
                        }
                    },
                    "varm": {
                        "membership": "attr",
                        "keys": {
                            "PCs": {
                                "membership": "item",
                            }
                        }
                    },
                    "uns": {
                        "membership": "attr",
                        "keys": {
                            "pca": {
                                "membership": "item",
                            }
                        }
                    },
                    "obsp": {
                        "membership": "attr",
                        "keys": {
                            "distances": {
                                "membership": "item",
                            },
                            "connectivities": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        gt_apis_dict["sc.tl.umap"]._products = BaseData()
        gt_apis_dict["sc.tl.umap"]._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obsm": {
                        "membership": "attr",
                        "keys": {
                            "X_umap": {
                                "membership": "item",
                            }
                        }
                    },
                    "uns": {
                        "membership": "attr",
                        "keys": {
                            "umap": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        gt_apis_dict["sc.tl.leiden"]._products = BaseData()
        gt_apis_dict["sc.tl.leiden"]._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obs": {
                        "membership": "attr",
                        "keys": {
                            "leiden": {
                                "membership": "item",
                            }
                        }
                    },
                    "uns": {
                        "membership": "attr",
                        "keys": {
                            "leiden": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        gt_apis_dict["sc.pl.umap"]._products = BaseData()
        gt_apis_dict["sc.pl.umap"]._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
            }
        )

        gt_apis_dict["root"]._products = BaseData()
        gt_apis_dict["root"]._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "X": {
                        "membership": "attr",
                    }
                }
            }
        )

        gt_deps_dict = {
            ("root", "sc.pp.neighbors"): BaseDependency(
                u_api_name="root",
                v_api_name="sc.pp.neighbors",
                args={
                    "data": "data"
                },
                arg_types={
                    "data": "object"
                },
            ),
            ("sc.pp.neighbors", "sc.tl.leiden"): BaseDependency(
                u_api_name="sc.pp.neighbors",
                v_api_name="sc.tl.leiden",
                args={
                    "data": "data"
                },
                arg_types={
                    "data": "object"
                },
            ),
            ("sc.pp.neighbors", "sc.tl.umap"): BaseDependency(
                u_api_name="sc.pp.neighbors",
                v_api_name="sc.tl.umap",
                args={
                    "data": "data"
                },
                arg_types={
                    "data": "object"
                },
            ),
            ("sc.tl.umap", "sc.pl.umap"): BaseDependency(
                u_api_name="sc.tl.umap",
                v_api_name="sc.pl.umap",
                args={
                    "data": "data"
                },
                arg_types={
                    "data": "object"
                },
            ),
            ("sc.tl.leiden", "sc.pl.umap"): BaseDependency(
                u_api_name="sc.tl.leiden",
                v_api_name="sc.pl.umap",
                args={
                    "color": "leiden"
                },
                arg_types={
                    "color": "str"
                },
            ),
        }

        gt_deps_dict[("root", "sc.pp.neighbors")].deps = BaseData()
        gt_deps_dict[("root", "sc.pp.neighbors")].deps.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "X": {
                        "membership": "attr",
                    }
                }
            }
        )

        gt_deps_dict[("sc.pp.neighbors", "sc.tl.umap")].deps = BaseData()
        gt_deps_dict[("sc.pp.neighbors", "sc.tl.umap")].deps.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "uns": {
                        "membership": "attr",
                        "keys": {
                            "neighbors": {
                                "membership": "item",
                            }
                        }
                    },
                    "obsp": {
                        "membership": "attr",
                        "keys": {
                            "connectivities": {
                                "membership": "item",
                            }
                        }
                    },
                }
            }
        )

        gt_deps_dict[("sc.pp.neighbors", "sc.tl.leiden")].deps = BaseData()
        gt_deps_dict[("sc.pp.neighbors", "sc.tl.leiden")].deps.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obsp": {
                        "membership": "attr",
                        "keys": {
                            "distances": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        gt_deps_dict[("sc.tl.umap", "sc.pl.umap")].deps = BaseData()
        gt_deps_dict[("sc.tl.umap", "sc.pl.umap")].deps.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obsm": {
                        "membership": "attr",
                        "keys": {
                            "X_umap": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        gt_deps_dict[("sc.tl.leiden", "sc.pl.umap")].deps = BaseData()
        gt_deps_dict[("sc.tl.leiden", "sc.pl.umap")].deps.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obs": {
                        "membership": "attr",
                        "keys": {
                            "leiden": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        assert graph.apis_dict == gt_apis_dict
        assert graph.deps_dict == gt_deps_dict

    def test_read_apis_from_graph_dict(self):
        mock_graph_dict = {
            "nodes": [
                {
                    "api": "sc.pp.neighbors",
                    "products": [
                        "data.obsm['X_pca']",
                        "data.varm['PCs']",
                        "data.uns['pca']",
                        "data.obsp['distances']",
                        "data.obsp['connectivities']"
                    ],
                    "id": "sc.pp.neighbors"
                },
                {
                    "api": "sc.pl.umap",
                    "products": [],
                    "id": "sc.pl.umap"
                },
                {   
                    "api": "root",
                    "products": [
                        "data.X"
                    ],
                    "id": "root"
                }
            ]
        }
        apis_dict = read_apis_from_graph_dict(mock_graph_dict, SCANPY_TOOLS_DICT)

        ScPpNeighbors = SCANPY_TOOLS_DICT["sc.pp.neighbors"]()
        ScPlUmap = SCANPY_TOOLS_DICT["sc.pl.umap"]()
        Root = SCANPY_TOOLS_DICT["root"]()

        ScPpNeighbors._api_name = "sc.pp.neighbors"
        ScPpNeighbors._products = BaseData()
        ScPpNeighbors._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obsm": {
                        "membership": "attr",
                        "keys": {
                            "X_pca": {
                                "membership": "item",
                            }
                        }
                    },
                    "varm": {
                        "membership": "attr",
                        "keys": {
                            "PCs": {
                                "membership": "item",
                            }
                        }
                    },
                    "uns": {
                        "membership": "attr",
                        "keys": {
                            "pca": {
                                "membership": "item",
                            }
                        }
                    },
                    "obsp": {
                        "membership": "attr",
                        "keys": {
                            "distances": {
                                "membership": "item",
                            },
                            "connectivities": {
                                "membership": "item",
                            }
                        }
                    }
                }
            }
        )

        ScPlUmap._api_name = "sc.pl.umap"
        ScPlUmap._products = BaseData()
        ScPlUmap._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
            }
        )

        Root._api_name = "root"
        Root._products = BaseData()
        Root._products.keys_info = BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "X": {
                        "membership": "attr",
                    }
                }
            }
        )

        assert apis_dict["sc.pp.neighbors"] == ScPpNeighbors
        assert apis_dict["sc.pl.umap"] == ScPlUmap
        assert apis_dict["root"] == Root
        
    def test_read_apis_from_graph_dict_1(self):
        """This is another test to mock more general cases."""
        mock_graph_dict = {
            "nodes": [
                {
                    "api": "api1",
                    "products": [ # complex nested cases
                        "data.obj1['obj2']['obj3']",
                        "data.obj1['obj2'].obj4",
                        "data['obj5'].obj6['obj7']",
                        "data['obj5'].obj6.obj8",
                    ],
                    "id": "api1"
                },
                {
                    "api": "api2",
                    "products": [ # replicated cases
                        "data.obj1['obj2']['obj3']",
                        "data.obj1['obj2']['obj3']",
                        "data.obj1['obj2']['obj3']",
                    ],
                    "id": "api2"
                },
                {
                    "api": "api3",
                    "products": [ # empty cases
                    ],
                    "id": "api3"
                }
            ]
        }

        apis_dict = read_apis_from_graph_dict(
            mock_graph_dict, 
            tools_dict={
                "api1": BaseAPI,
                "api2": BaseAPI,
                "api3": BaseAPI,
            }
        )

        assert apis_dict["api1"]._products.keys_info == BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obj1": {
                        "membership": "attr",
                        "keys": {
                            "obj2": {
                                "membership": "item",
                                "keys": {
                                    "obj3": {
                                        "membership": "item",
                                    },
                                    "obj4": {
                                        "membership": "attr",
                                    },
                                },
                            },
                        },
                    },
                    "obj5": {
                        "membership": "item",
                        "keys": {
                            "obj6": {
                                "membership": "attr",
                                "keys": {
                                    "obj7": {
                                        "membership": "item",
                                    },
                                    "obj8": {
                                        "membership": "attr",
                                    },
                                },
                            },
                        },
                    },
                }
            }
        )

        assert apis_dict["api2"]._products.keys_info == BaseKeysInfo.model_validate(
            {
                "membership": "self",
                "keys": {
                    "obj1": {
                        "membership": "attr",
                        "keys": {
                            "obj2": {
                                "membership": "item",
                                "keys": {
                                    "obj3": {
                                        "membership": "item",
                                    },
                                },
                            },
                        },
                    },
                },
            }
        )

        assert apis_dict["api3"]._products.keys_info == BaseKeysInfo.model_validate(
            {
                "membership": "self",
            }
        )
        
    def test_read_deps_from_graph_dict(self):
        pass

    def test_retrieve_products(self):
        test_data = TestData()
        data = BaseData(
            data=test_data,
            keys_info=test_data.keys_info
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
            data=BaseTestData(),
            keys_info=dep_keys_info
        )

        dependency = BaseDependency(
            u_api_name="src",
            v_api_name="dst",
            args={},
            arg_types={},
            deps=dep_data,
        )

        result = retrieve_products(src_api, dependency)

        assert result.deps.data.obj1 == src_api._products.data.obj1
        assert result.deps.data.obj2["key1"] == src_api._products.data.obj2["key1"]
        assert "key2" not in result.deps.data.obj2.keys()
        assert result.deps.data["key3"].key4 == src_api._products.data["key3"].key4

    def test_aggregate_deps(self):
        test_data = TestData()
        test_data_src1 = TestDataSrc1()
        test_data_src2 = TestDataSrc2()
        test_data_gt = TestDataGT()
        dep1 = BaseDependency(
            u_api_name="src1",
            v_api_name="dst",
            args={},
            arg_types={},
            deps=BaseData(
                data=test_data_src1,
                keys_info=test_data_src1.keys_info
            ),
        )

        dep2 = BaseDependency(
            u_api_name="src2",
            v_api_name="dst",
            args={},
            arg_types={},
            deps=BaseData(
                data=test_data_src2,
                keys_info=test_data_src2.keys_info
            ),
        )

        
        data = BaseData(
            data=test_data,
            keys_info=test_data.keys_info
        )

        dst_api = BaseAPI()
        dst_api._deps = data

        dst_api = aggregate_deps([dep1, dep2], dst_api)

        gt_data = BaseData(
            data=test_data_gt,
            keys_info=test_data_gt.keys_info
        )

        assert dst_api._deps.data == gt_data.data
        assert dst_api._deps.keys_info == gt_data.keys_info

if __name__ == "__main__":
    """__main__ is created for easy debugging without running the whole test suite."""
    test_case = TestDepGraph()
    test_case.test_read_apis_from_graph_dict()
    test_case.test_read_apis_from_graph_dict_1()
    test_case.test_retrieve_products()
    test_case.test_construct_dep_graph()
    test_case.test_aggregate_deps()
    unittest.main()
    # pass