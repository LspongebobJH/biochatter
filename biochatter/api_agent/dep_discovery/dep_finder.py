import logging
logger = logging.getLogger(__name__)

from typing import Literal
from importlib import import_module
import networkx as nx
import libcst as cst

from .utils import set_state, build_API_calling, Knockoff, CodeHistory
from biochatter.api_agent.base.utils import run_codes

class DependencyFinder:
    MAX_CNT = 3
    SEP_LINE = "-" * 40
    # Jiahang: remove these constants.
    _START_EACH_CHAIN = 0
    _FINISH_OBJ_DEP_EACH_CHAIN = _START_EACH_CHAIN + 1
    _FINISH_API_DEP_EACH_CHAIN = _FINISH_OBJ_DEP_EACH_CHAIN + 1
    _FINISH_API_DEP_ALL_CHAINS = _FINISH_API_DEP_EACH_CHAIN + 1
    _FINISH_API_DEP_GRAPH = _FINISH_API_DEP_ALL_CHAINS + 1

    def __init__(self, 
                 api_set_name: str,
                 api_sel: str | list[int], 
                 state_name: dict,
                 api_dict: dict,
                 repr_type: Literal['nest_dict', 'code_line'] = 'code_line',
                 debug: bool = False):
        """
        Discovery dependency graph among API chains given multiple API chains represented in code lines, where
        each APi chain is seens as a use case of the package.

        Parameters
        ----------
        api_set_name : str
            The API set path relative to 'resource' directory, such as 'scanpy.SCANPY_PL_MULTI_LINE_SIMPLE_ARG'.
        api_sel : str | int
            The slice of API selections. If 'all' then select all API. If a list of integers, then select the API indexed by `api_sel`.
        state_name : str
            The state name of the code execution env. Available states are ['scanpy', 'squidpy'].
        repr_type : Literal['nest_dict', 'code_line'], default 'code_line'
            The representation type of the API. If 'nest_dict', then the API is represented as a nested dictionary.
            If 'code_line', then the API is represented as a code line. For now, only 'code_line' is supported.
        debug : bool, default False
            Whether to run in debug mode. If True, then only run the first 3 API chains within the selected API subset.
        """
        # Jiahang(severe): this is a temporary solution to import the API chains.
        self.api_dict = api_dict
        self.api_chains = import_module(f'biochatter.api_agent.dep_discovery.SCANPY_PL_MULTI_LINE_SIMPLE_ARG').DATA
        self.api_meta = import_module(f'biochatter.api_agent.dep_discovery.SCANPY_PL_MULTI_LINE_SIMPLE_ARG').META
        self.api_sel = range(len(self.api_chains)) if api_sel == 'all' else api_sel
        self.repr_type = repr_type
        assert self.repr_type == 'code_line', "Only code_line representation is supported for now."
        self.state_name = state_name
        self.debug = debug

        self.code_history = CodeHistory()
        self.knockoff_runner = Knockoff()

        self.run_level = DependencyFinder._START_EACH_CHAIN # this is used to easy track execution flow.

    def clear(self):
        self.code_history.clear()
        self.knockoff_runner.clear()
        self.run_level = DependencyFinder._START_EACH_CHAIN


    def find_obj_dep_for_each_chain(self, api_chain: list):
        assert self.run_level < DependencyFinder._FINISH_OBJ_DEP_EACH_CHAIN, \
            f"Execution flow should not finish object-level dependency discovery for each chain." \
            f"Please check the run level {self.run_level}."
    
        num_api = len(api_chain)
        res_api_chain = []
        all_dep_prod = []
        state = set_state(self.state_name)

        res_api_chain.append(api_chain[0])

        for dep_idx in range(num_api-1):
            dep_api = api_chain[dep_idx] # the dependent api
            target_api = api_chain[dep_idx + 1] # the next api

            dep_api_call = build_API_calling(dep_api)
            self.code_history.clear()
            self.code_history.add_code(dep_api_call)
            
            output, error = run_codes(self.code_history.code_history, state=state)
            if isinstance(error, Exception):
                logger.info(f"Program failed in initial calling: {output}\n"
                            f"Codes:\n{self.code_history.code_history}")
                return None, False
            
            dep_prod = dep_api['products']
            all_dep_prod.extend(dep_prod)

            """Key: discover unitary dependency of each target API on products produced by the all dependent(previous) API"""
            self.knockoff_runner.preprocess(all_dep_prod, target_api, self.state_name, state)
            nec_prod = self.knockoff_runner.run_unitary()

            target_api['dependencies_obj'] = nec_prod
            res_api_chain.append(target_api)

        self.run_level = DependencyFinder._FINISH_OBJ_DEP_EACH_CHAIN
        return res_api_chain, True
    
    def find_api_dep_for_each_chain(self, res_api_chain: list):
        assert self.run_level < DependencyFinder._FINISH_API_DEP_EACH_CHAIN, \
            f"Execution flow should not finish API-level dependency discovery for each chain." \
            f"Please check the run level {self.run_level}."
        num_api = len(res_api_chain)
        _res_api_chain = res_api_chain.copy()

        for cur_idx in range(1, num_api):
            cur_api: dict = res_api_chain[cur_idx]
            cur_dependencies = cur_api.pop('dependencies_obj')
            cur_dependencies_api = {}
            for dep_idx in range(cur_idx):
                dep_api = res_api_chain[dep_idx]
                dep_products = dep_api['products']
                consumed_prod = list(set(cur_dependencies).intersection(set(dep_products)))
                if consumed_prod: # has API dependency
                    cur_dependencies_api[dep_api['api']] = consumed_prod
            cur_api['dependencies'] = cur_dependencies_api
            _res_api_chain[cur_idx] = cur_api # this line is indeed unnecessary, but still used here.

        print(DependencyFinder.SEP_LINE)

        res_api_chain = _res_api_chain
        self.run_level = DependencyFinder._FINISH_API_DEP_EACH_CHAIN
        return res_api_chain
    
    def _parse_codes(self, codes: str) -> dict:
        code_lines = [line.strip() for line in codes.split('\n')]
        code_lines = [line for line in code_lines if line]

        results = []

        # Find the function call
        class FunctionCallVisitor(cst.CSTVisitor):
            def __init__(self):
                self.func_name = None
                self.args = {}
                self.arg_types = {}
            
            def visit_Call(self, node: cst.Call):
                # Get function name
                if isinstance(node.func, cst.Name):
                    self.func_name = node.func.value
                elif isinstance(node.func, cst.Attribute):
                    parts = []
                    current = node.func
                    while isinstance(current, cst.Attribute):
                        parts.append(current.attr.value)
                        current = current.value
                    if isinstance(current, cst.Name):
                        parts.append(current.value)
                    self.func_name = ".".join(reversed(parts))
                
                # Get arguments
                for arg in node.args:
                    if arg.keyword:
                        # Determine argument type
                        if isinstance(arg.value, cst.Name):
                            self.arg_types[arg.keyword.value] = "object"
                        elif isinstance(arg.value, cst.Integer):
                            self.arg_types[arg.keyword.value] = "int"
                        elif isinstance(arg.value, cst.Float):
                            self.arg_types[arg.keyword.value] = "float"
                        elif isinstance(arg.value, cst.SimpleString):
                            self.arg_types[arg.keyword.value] = "str"
                        elif isinstance(arg.value, cst.List):
                            self.arg_types[arg.keyword.value] = "list"
                        elif isinstance(arg.value, cst.Tuple):
                            self.arg_types[arg.keyword.value] = "tuple"
                        elif isinstance(arg.value, cst.Dict):
                            self.arg_types[arg.keyword.value] = "dict"
                        elif isinstance(arg.value, cst.Set):
                            self.arg_types[arg.keyword.value] = "set"
                        elif arg.value.value in ['True', 'False']:
                            self.arg_types[arg.keyword.value] = "bool"
                        else:
                            self.arg_types[arg.keyword.value] = "object"
                        
                        # when str, arg.value.value is tricky, it would be a string of string.
                        self.args[arg.keyword.value] = arg.value.value

        for code_line in code_lines:
            # Parse the code line using libcst
            module = cst.parse_module(code_line)
            visitor = FunctionCallVisitor()
            module.visit(visitor)
            
            if visitor.func_name:
                results.append({
                    'api': visitor.func_name,
                    'args': visitor.args,
                    'arg_types': visitor.arg_types
                })

        return results

    def _add_products(self, parsed_api_chain: list) -> list:
        # import products from pre-compiled API data models
        for idx in range(len(parsed_api_chain)):
            parsed_api_chain[idx]['products'] = self.api_dict[parsed_api_chain[idx]['api']]['products']
            
        return parsed_api_chain

    def find_obj_dep_for_all_chains(self):
        assert self.run_level == DependencyFinder._START_EACH_CHAIN, \
            f"Execution flow should be at the start status. Please check the run level {self.run_level}."
        cnt = 0
        res_api_chains = []
        
        for api_chain in [self.api_chains[i] for i in self.api_sel]:
            self.run_level = DependencyFinder._START_EACH_CHAIN
            codes = api_chain['codes']
            parsed_api_chain = self._parse_codes(codes)
            parsed_api_chain = self._add_products(parsed_api_chain)
            res_api_chain, success = self.find_obj_dep_for_each_chain(parsed_api_chain)
            if not success:
                continue    
            res_api_chain = self.find_api_dep_for_each_chain(res_api_chain)
            res_api_chains.append(res_api_chain)
            cnt += 1

            if self.debug and cnt > DependencyFinder.MAX_CNT:
                break
        self.run_level = DependencyFinder._FINISH_API_DEP_ALL_CHAINS
        return res_api_chains
    
    def construct_dep_graph(self): # the entry point of the whole workflow
        res_api_chains = self.find_obj_dep_for_all_chains()
        G = nx.DiGraph()
        for api_chain in res_api_chains:
            for api in api_chain:
                _api: dict = api.copy()
                dependencies: dict = _api.pop('dependencies') if 'dependencies' in _api.keys() else {}
                G.add_node(_api['api'], **_api) # dependencies should not be added
                for dep_api_name, obj_dep in dependencies.items():
                    G.add_edge(dep_api_name, _api['api'], dependencies=obj_dep)

        # create a root node
        G.add_node('root')
        in_degrees = G.in_degree()
        init_nodes = [
            node for node, degree in in_degrees if degree == 0 and node != 'root'
        ]
        G.add_edges_from([('root', node) for node in init_nodes])

        self.run_level = DependencyFinder._FINISH_API_DEP_GRAPH
        return G, res_api_chains
    
    def __call__(self, *args, **kwds):
        return self.construct_dep_graph()
        


            
        
    

    

