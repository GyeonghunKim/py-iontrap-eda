from typing import Dict
from .layer.base_layer import LayerName
from shapely.geometry import Polygon

class Electrode:
    def __init__(self, name: str, geometry: Dict[LayerName, Polygon]):
        self.name = name
        self.geometry = geometry

    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"Electrode(name={self.name}, geometry={self.geometry})"
    
    