from typing import Tuple, Optional
import shapely

from .route_method import RouteMethod
from .port import Port
from .layer.base_layer import LayerName

class BondPad:
    def __init__(self, 
                 name: str, 
                 center: Tuple[float, float], 
                 width: float, 
                 height: float,
                 port_orientation: int,
                 route_method: RouteMethod=RouteMethod.VIA,
                 via_pad_center: Optional[Tuple[float, float]]=None,
                 via_pad_width: Optional[float]=None,
                 via_pad_height: Optional[float]=None,
                 via_pad_layer_down_to: Optional[LayerName]=None,
                 ):
        self.name = name
        self.center = center
        self.width = width
        self.height = height
        self.port_orientation = port_orientation
        self.geometry = shapely.Polygon(
            [(center[0] - width/2, center[1] - height/2),
             (center[0] - width/2, center[1] + height/2),
             (center[0] + width/2, center[1] + height/2),
             (center[0] + width/2, center[1] - height/2)]
        )
        self.route_method = route_method
        if route_method == RouteMethod.VIA:
            if via_pad_center is None or via_pad_width is None or via_pad_height is None:
                raise ValueError("via_pad_center, via_pad_width, and via_pad_height must be provided for via route method")
            self.via_pad = shapely.Polygon(
                [(via_pad_center[0] - via_pad_width/2, via_pad_center[1] - via_pad_height/2),
                 (via_pad_center[0] - via_pad_width/2, via_pad_center[1] + via_pad_height/2),
                 (via_pad_center[0] + via_pad_width/2, via_pad_center[1] + via_pad_height/2),
                 (via_pad_center[0] + via_pad_width/2, via_pad_center[1] - via_pad_height/2)]
            )
            self.via_pad_layer_down_to = via_pad_layer_down_to
        port_geometry = shapely.LineString()
        if port_orientation == 0:
            port_geometry = shapely.LineString([(center[0] + width/2, center[1] - height/2), (center[0] + width/2, center[1] + height/2)])
        elif port_orientation == 180:
            port_geometry = shapely.LineString([(center[0] - width/2, center[1] - height/2), (center[0] - width/2, center[1] + height/2)])
        elif port_orientation == 90:
            port_geometry = shapely.LineString([(center[0] - width/2, center[1] + height/2), (center[0] + width/2, center[1] + height/2)])
        elif port_orientation == -90:
            port_geometry = shapely.LineString([(center[0] - width/2, center[1] - height/2), (center[0] + width/2, center[1] - height/2)])
        else:
            raise ValueError(f"Invalid port orientation: {port_orientation}")
        self.port = Port(name, port_geometry, port_orientation)
        
    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"BondPad(name={self.name}, center={self.center}, width={self.width}, height={self.height}, port_orientation={self.port_orientation})"

    def __eq__(self, other):
        return self.name == other.name and self.center == other.center and self.width == other.width and self.height == other.height and self.port_orientation == other.port_orientation

    def __hash__(self):
        return hash((self.name, self.center, self.width, self.height, self.port_orientation))