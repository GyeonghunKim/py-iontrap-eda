from typing import List, Optional

import numpy as np
import shapely


class Zone:
    def __init__(self, name: str, dc_widths: np.ndarray, negatives: Optional[List[shapely.Polygon]]=None):
        self.name = name
        self.dc_widths = dc_widths
        self.negatives: List[shapely.Polygon] = negatives if negatives is not None else []
        self.center_z: float = 0
        
    def translate(self, d_center_z: float):
        self.center_z += d_center_z
        
    @property
    def width(self):
        return np.sum(self.dc_widths)


class Spacing:
    def __init__(self, width: Optional[float] = None, n_electrodes: int = 1, negatives: Optional[List[shapely.Polygon]]=None):
        self.width = width
        self.n_electrodes = n_electrodes
        self.center_z: float = 0
        self.negatives: List[shapely.Polygon] = negatives if negatives is not None else []

    def set_width(self, width: float):
        self.width = width
        
    def translate(self, d_center_z: float):
        if self.center_z is None:
            raise ValueError("Cannot translate spacing with undefined width")
        self.center_z += d_center_z
        
    @property
    def dc_widths(self):
        return np.array([self.width / self.n_electrodes] * self.n_electrodes)
