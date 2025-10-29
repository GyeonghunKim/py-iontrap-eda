from typing import List, Dict
from .layer.base_layer import BaseLayer, LayerName
from .thickness_map import ThicknessMap
from .rail.base_rail import BaseRail
class Trap:
    def __init__(self, 
                 name: str, 
                 description: str, 
                 width: float,
                 height: float,
                 dicing_margin: float,
                 dicing_channel_layer: LayerName,
                 layers: Dict[LayerName, BaseLayer], 
                 thickness_map: ThicknessMap):
        self.name = name
        self.description = description
        self.width = width
        self.height = height
        self.dicing_margin = dicing_margin
        self.layers = layers
        self.thickness_map = thickness_map  
        
        

    def add_rail(self, rail: BaseRail):
        outer_layer_name = self.thickness_map.layer_names[0]
        outer_layer = self.layers[outer_layer_name]

        # for electrode in rail.electrodes:
            