from typing import Optional, Tuple
import shapely
import numpy as np

from .port import Port

from .route_method import RouteMethod
    
class Electrode:
    def __init__(self, name: str, geometry: shapely.Polygon, route_method: RouteMethod=RouteMethod.VIA):
        self.name = name
        self.geometry = geometry
        self.port: Optional[Port] = None
        self.via_pad: Optional[shapely.Polygon] = None
        self.route_method = route_method
        self.via_pad_layer_down_to: Optional[LayerName] = None
    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"Electrode(name={self.name}, geometry={self.geometry})"

    def change_via_port_orientation(self, orientation: int):
        if self.port is None:
            raise ValueError(f"Port is not set for electrode {self.name}")
        if self.port.orientation == orientation:
            return
        min_x, min_y, max_x, max_y = self.via_pad.bounds    
    
        if orientation == 0:
            self.port = Port(self.name, shapely.LineString([(max_x, min_y), (max_x, max_y)]), orientation)
        elif orientation == 180:
            self.port = Port(self.name, shapely.LineString([(min_x, min_y), (min_x, max_y)]), orientation)
        elif orientation == 90:
            self.port = Port(self.name, shapely.LineString([(min_x, max_y), (max_x, max_y)]), orientation)
        elif orientation == -90:
            self.port = Port(self.name, shapely.LineString([(min_x, min_y), (max_x, min_y)]), orientation)
        else:
            raise ValueError(f"Invalid orientation: {orientation}")