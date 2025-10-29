from .base_layer import BaseLayer, LayerName
from typing import Dict
from .layer_electrode import LayerElectrode

class ViaLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str, thickness: float, electrodes: Dict[str, LayerElectrode]):
        super().__init__(name, long_name, description, thickness, electrodes)