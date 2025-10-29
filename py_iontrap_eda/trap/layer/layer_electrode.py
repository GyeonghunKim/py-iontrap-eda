from typing import Dict, List
from shapely.geometry import Polygon

class LayerElectrode:
    def __init__(self, name: str, geometry: List[Polygon]):
        self.name = name
        self.geometry = geometry

    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"LayerElectrode(name={self.name}, geometry={self.geometry})"
    
    