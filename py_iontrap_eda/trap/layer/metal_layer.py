from .base_layer import BaseLayer, LayerName
from typing import Dict, List, Union
from ..electrode import Electrode
import shapely
class MetalLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str):
        super().__init__(name, long_name, description)
        self.total_negatives = []
        self.electrodes: Dict[str, Electrode] = {}  
    def add_total_negatives(self, total_negatives: List[shapely.Polygon]):
        self.total_negatives = self.total_negatives + total_negatives
        
    def add_electrode(self, electrode: Electrode):
        self.electrodes[electrode.name] = electrode
        