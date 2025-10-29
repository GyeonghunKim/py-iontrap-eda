from typing import NewType

LayerName = NewType('LayerName', int)

class BaseLayer:
    def __init__(self, name: LayerName, long_name: str, description: str, thickness: float, electrodes: Dict[str, LayerElectrode]):
        self.name = name
        self.long_name = long_name
        self.description = description
        self.thickness = thickness
        self.electrodes = electrodes

    def __str__(self):
        return f"{self.name}: {self.long_name} ({self.description})"

    def __repr__(self):
        return f"BaseLayer(name={self.name}, long_name={self.long_name}, description={self.description}, thickness={self.thickness})"
    
    