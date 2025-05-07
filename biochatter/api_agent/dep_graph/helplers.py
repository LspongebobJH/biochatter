from pydantic import BaseModel

class TopoSortLayers(list):
    def __init__(self, layers: list[list[BaseModel]]):
        self.layers = layers
        self.index = 0

    def update_layer(self, idx: int, new_layer: list[BaseModel]):
        self.layers[idx] = new_layer
    
    def __getitem__(self, idx: int) -> list[BaseModel]:
        return self.layers[idx]
    
    @classmethod
    def retrieve_apis_if_in_layer(cls, apis: list[BaseModel], layer: list[BaseModel]) -> list[BaseModel]:
        """Retrieve apis in given apis if in the given layer."""
        return [api for api in apis if api in layer]