from .base_layer import BaseLayer, LayerName
from typing import Dict
from ..electrode import Electrode

class ViaLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str):
        super().__init__(name, long_name, description)