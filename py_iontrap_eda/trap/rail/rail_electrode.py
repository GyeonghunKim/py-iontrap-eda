from typing import Optional, Tuple
from enum import Enum, auto
import shapely
import numpy as np

from .port import Port

class RouteMethod(Enum):
    VIA = auto()
    PORT = auto()
    
class RailElectrode:
    def __init__(self, name: str, geometry: shapely.Polygon, route_method: RouteMethod=RouteMethod.VIA):
        self.name = name
        self.geometry = geometry
        self.port: Optional[Port] = None
        self.via_pad: Optional[shapely.Polygon] = None
        self.route_method = route_method
        
    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"RailElectrode(name={self.name}, geometry={self.geometry})"
