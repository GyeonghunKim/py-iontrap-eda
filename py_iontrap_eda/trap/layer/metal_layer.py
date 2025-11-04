from .base_layer import BaseLayer, LayerName
from typing import Dict, List, Union
from ..electrode import Electrode
import shapely
class MetalLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str):
        super().__init__(name, long_name, description)
        self.total_negatives = []
        self.electrodes: Dict[str, Electrode] = {}  
        self.geometries: List[shapely.Polygon] = []
        self.via_pads: List[shapely.Polygon] = []
        
    def add_total_negatives(self, total_negatives: List[shapely.Polygon]):
        self.total_negatives = self.total_negatives + total_negatives
        
    def add_electrode(self, electrode: Electrode):
        self.electrodes[electrode.name] = electrode
        
    def add_via_pad(self, via_pad: shapely.Polygon):
        self.via_pads.append(via_pad)
        
    def add_via_pads(self, via_pads: List[shapely.Polygon]):
        self.via_pads = self.via_pads + via_pads
        
    def __str__(self):
        return f"MetalLayer(name={self.name}, electrodes={self.electrodes})"
    
    def __repr__(self):
        return f"MetalLayer(name={self.name}, electrodes={self.electrodes}, via_pads={self.via_pads})"
        