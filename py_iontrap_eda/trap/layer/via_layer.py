from .base_layer import BaseLayer, LayerName
from typing import Dict
from ..electrode import Electrode

class ViaLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str):
        super().__init__(name, long_name, description)
        self.via_areas: List[shapely.Polygon] = []
        
    def set_via_dimensions(self, width: float, height: float, spacing: float):
        self.width = width
        self.height = height
        self.spacing = spacing
        
    def add_via_area(self, via_area: shapely.Polygon):
        self.via_areas.append(via_area)
        
    def add_via_areas(self, via_areas: List[shapely.Polygon]):
        self.via_areas = self.via_areas + via_areas
        