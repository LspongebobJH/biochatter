import logging

logging.basicConfig(
    level=logging.INFO,
    # filename=f'exp/dep_graph/logs/scanpy_api_dep/{api_type}_{date}.log',
    format='%(asctime)s - %(name)s - %(levelname)s \n%(message)s\n',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

from .dep_finder import DependencyFinder
import datetime
import json
import networkx as nx
import os
import pickle
import argparse

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--debug', action='store_true', default=False, help='Debug mode')
    # args = parser.parse_args()
    
    """Initialization for whole experiment"""
    api_set_name = 'scanpy.SCANPY_PL_MULTI_LINE_SIMPLE_ARG'
    debug = True
    api_sel = [0]
    date = datetime.datetime.now().strftime("%d-%m-%y")
    result_dir = f'exp/dep_graph/results/scanpy_api_dep/{api_set_name}_{date}/'

    dep_finder = DependencyFinder(
        api_set_name=api_set_name,
        api_sel=api_sel,
        state_name='scanpy',
        debug=debug
    )

    G, res_APIs = dep_finder()

    """Saving results"""
    if not debug:
        os.makedirs(os.path.dirname(result_dir), exist_ok=True)
        data_path = os.path.join(result_dir, 'data.json')
        graph_path = os.path.join(result_dir, 'graph.pkl')
        graph_json_path = os.path.join(result_dir, 'graph.json')

        with open(data_path, 'w') as f:
            json.dump(res_APIs, f, indent=4)
        logger.info(f"Results saved to {data_path}")

        with open(graph_path, 'wb') as f:
            pickle.dump(G, f)
        logger.info(f"Graph saved to {graph_path}")

        with open(graph_json_path, 'w') as f:
            json.dump(nx.node_link_data(G, edges='edges'), f, indent=4)
        logger.info(f"Graph json saved to {graph_json_path}")
    
    pass