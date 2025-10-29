from typing import NewType, Tuple, Dict
from .layer_electrode import LayerElectrode

LayerName = NewType('LayerName', Tuple[int, int])

class BaseLayer:
    def __init__(self, name: LayerName, long_name: str, description: str):
        self.name = name
        self.long_name = long_name
        self.description = description
        self.electrodes: List[LayerElectrode] = []

    def __str__(self):
        return f"{self.name}: {self.long_name} ({self.description})"

    def __repr__(self):
        return f"BaseLayer(name={self.name}, long_name={self.long_name}, description={self.description}, thickness={self.thickness})"
    
    