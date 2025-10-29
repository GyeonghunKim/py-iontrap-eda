from .base_layer import BaseLayer, LayerName
from typing import Dict
from .layer_electrode import LayerElectrode

class MetalLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str):
        super().__init__(name, long_name, description)