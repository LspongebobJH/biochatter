from unittest.mock import MagicMock

from biochatter.api_agent.base.agent_abc import BaseAPI, BaseDependency, BaseData, BaseKeysInfo
from biochatter.api_agent.base.utils import retrieve_products, read_apis_from_graph_dict, read_deps_from_graph_dict
from biochatter.api_agent.python.scanpy import SCANPY_TOOLS_DICT
from biochatter.api_agent.dep_graph.graph import DependencyGraph
def test_construct_dep_graph():
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

    pass



def test_read_apis_from_graph_dict():
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
    

def test_read_apis_from_graph_dict_1():
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
    
def test_read_deps_from_graph_dict():
    pass

def test_retrieve_products():
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
    test_read_apis_from_graph_dict()
    test_read_apis_from_graph_dict_1()
    test_retrieve_products()
    test_construct_dep_graph()
    pass