from dataclasses import dataclass
from typing import List, Optional
import pickle



@dataclass
class RadialDesign:
    inner_dc_width: float
    rf_width: float

    def save(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
    @classmethod    
    def load(cls, filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)

@dataclass
class AxialDesign:
    offset_from_center: float
    n_participation_electrode_pairs: int
    widths: List[float]
    
@dataclass
class TweakerDesign:
    offset_from_center: float
    width: float
    
    
@dataclass
class TrapDesign:
    name: str
    radial_design: RadialDesign
    axial_design: AxialDesign
    tweaker_design: Optional[TweakerDesign] = None
    