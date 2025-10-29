from typing import List
from ..layer import Layer
from .thickness_map import ThicknessMap
from .units import Units
from .rail import Rail
class Trap:
    def __init__(self, name: str, description: str, layers: Dict[LayerName, Layer], thickness_map: ThicknessMap):
        self.name = name
        self.description = description
        self.layers = layers
        self.thickness_map = thickness_map  

    def add_rail(self, rail: Rail):
        outer_layer_name = self.thickness_map.layer_names[0]
        outer_layer = self.layers[outer_layer_name]

        for electrode in rail.electrodes:
            