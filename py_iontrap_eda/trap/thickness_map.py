from typing import List
from .layer.base_layer import LayerName

class ThicknessMap:
    def __init__(self, layer_names: List[LayerName], layer_thicknesses: List[float]):
        self.layer_names = layer_names
        self.layer_thicknesses = layer_thicknesses

    def __str__(self):
        return "".join([f"{layer_name} -> {layer_thickness}\n" for layer_name, layer_thickness in zip(self.layer_names, self.layer_thicknesses)])

    def __repr__(self):
        return f"ThicknessMap(layer_names={self.layer_names}, layer_thicknesses={self.layer_thicknesses})"
        
    def get_thickness(self, layer_name: LayerName) -> float:
        return self.layer_thicknesses[self.layer_names.index(layer_name)]