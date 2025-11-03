from typing import List, Union, Dict
from abc import abstractmethod, ABC
from dataclasses import dataclass, field
import pickle

import numpy as np
import shapely
import matplotlib.pyplot as plt
import gdsfactory as gf

from ..electrode import Electrode
from ..route_method import RouteMethod
from .zone import Zone, Spacing
from ..port import Port
from ...utils import Units

@dataclass
class BaseRailParameters:
    dc_length: float
    total_width: float
    zones: List[Union[Zone, Spacing]]
    z_padding: float
    via_pad_width: float

class BaseRail(ABC):
    def __init__(self, name: str, parameters: BaseRailParameters):
        self.name = name
        self.parameters = parameters
        self.electrodes: List[Electrode] = []
        self._assign_width_blank_spaces()
        self.dc_widths = np.hstack([zones.dc_widths for zones in self.parameters.zones])
        self.dc_z_points = np.hstack(
            [
                [-self.parameters.total_width / 2],
                np.cumsum(self.dc_widths) - self.parameters.total_width / 2,
            ]
        )
        self.center_z = parameters.total_width / 2
        self._translate_zones()
        self.total_negatives = self._get_total_negatives()
    
    def get_electrode(self, name: str):
        for electrode in self.electrodes:
            if electrode.name == name:
                return electrode
        return None
    
    def save(self, filename: str):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
    @classmethod
    def load(cls, filename: str):
        with open(filename, 'rb') as f:
            return pickle.load(f)

    def set_route_method(self, method_dict: Dict[str, RouteMethod]):
        for electrode in self.electrodes:
            if electrode.name in method_dict:
                electrode.route_method = method_dict[electrode.name]
            
    def _assign_width_blank_spaces(self):
        blank_space_electrodes = 0
        occupied_width = 0

        for zone in self.parameters.zones:
            if isinstance(zone, Zone):
                occupied_width += np.sum(zone.dc_widths)
            elif isinstance(zone, Spacing):
                if zone.width is None:
                    blank_space_electrodes += zone.n_electrodes
                else:
                    occupied_width += zone.width
            else:
                raise TypeError(
                    f"Zone must be of type Zone or Spacing, got {type(zone)}"
                )

        if blank_space_electrodes == 0:
            return

        blank_electrode_width = (
            self.parameters.total_width - occupied_width
        ) / blank_space_electrodes

        if blank_electrode_width < 0:
            raise ValueError(
                f"Total width of electrodes ({occupied_width:.2f}) exceeds rail width ({self.parameters.total_width:.2f})"
            )

        for zone in self.parameters.zones:
            if isinstance(zone, Spacing):
                if zone.width is None:
                    zone.set_width(zone.n_electrodes * blank_electrode_width)
                    
    def _translate_zones(self):
        for i, zone in enumerate(self.parameters.zones):
            if i == 0:
                zone.translate(-self.parameters.total_width / 2 + zone.width / 2)
            else:
                prev_zones = self.parameters.zones[:i]
                prev_width = sum(
                    np.sum(z.dc_widths) for z in prev_zones
                )
                zone.translate(-self.parameters.total_width / 2 + prev_width + zone.width / 2)

    def get_bbox(self):
        return shapely.box(
            *shapely.MultiPolygon(
                [electrode.geometry for electrode in self.electrodes]
            ).bounds
        )
        
    def translate(self, x: float, y: float):
        self.center_z += x
        for electrode in self.electrodes:
            electrode.geometry = shapely.affinity.translate(electrode.geometry, x, y)
        self.total_negatives = shapely.affinity.translate(
            self.total_negatives, x, y
            )
        self.generate_ports()
        
    def generate_ports(self):
        bbox = self.get_bbox()
        bbox_x_min, bbox_y_min, bbox_x_max, bbox_y_max = bbox.bounds
        for electrode in self.electrodes:
            intersection = bbox.exterior.intersection(electrode.geometry.exterior)
            if intersection.geom_type == 'LineString':
                x_min, y_min, x_max, y_max = intersection.bounds
                if x_min == x_max: # either 0 or 180
                    if np.isclose(x_min, bbox_x_min):
                        orientation = 180 
                    elif np.isclose(x_min, bbox_x_max):
                        orientation = 0
                    else:
                        raise ValueError(f"Expected LineString intersection for electrode {electrode.name}, got {intersection.geom_type}")
                elif y_min == y_max: # either 90 or -90
                    if np.isclose(y_min, bbox_y_min):
                        orientation = 90
                    elif np.isclose(y_min, bbox_y_max):
                        orientation = -90
                    else:
                        raise ValueError(f"Expected LineString intersection for electrode {electrode.name}, got {intersection.geom_type}")
                else:
                    raise ValueError(f"Expected LineString intersection for electrode {electrode.name}, got {intersection.geom_type}")
                electrode.port = Port(electrode.name, intersection, orientation)
            else:
                raise ValueError(f"Exxpected LineString intersection for electrode {electrode.name}, got {intersection.geom_type}")
        
    def generate_via_pads(self):
        for electrode in self.electrodes:
            electrode.via_pad = shapely.intersection(
                    electrode.geometry, 
                    electrode.port.geometry.buffer(self.parameters.via_pad_width, cap_style='square')
                )
     

    def _get_total_negatives(self):
        negatives = []
        for zone in self.parameters.zones:
            for negative in zone.negatives:
                negatives.append(shapely.affinity.translate(negative, zone.center_z))
        total_negative_geom = shapely.union_all(negatives)
        return total_negative_geom
    
    def show(self, show_via_pads: bool=True): # , fig=None, ax=None, show_bbox: bool=True):
        c = gf.Component()
        n = gf.Component()
        for i, negative in enumerate(self.total_negatives.geoms):
            n.add_polygon(np.array(negative.exterior.xy).T / Units.um, layer=(1, 0))
        
        for i, electrode in enumerate(self.electrodes):
            c_temp = gf.Component()
            c_temp.add_polygon(np.array(electrode.geometry.exterior.xy).T / Units.um, layer=(1, 0))
            c_temp = gf.boolean(c_temp, n, 'A-B', layer=(1, 0))
            c_temp.add_label(f'{electrode.name}', position=c_temp.dcenter, layer=(1, 0))
            c_temp.add_port(name=electrode.port.name, center=np.array(electrode.port.center) / Units.um, width=electrode.port.width / Units.um, orientation=electrode.port.orientation, layer=(1, 0))
            c << c_temp

        if show_via_pads:
            for electrode in self.electrodes:
                if electrode.via_pad is not None:
                    c.add_polygon(np.array(electrode.via_pad.exterior.xy).T / Units.um, layer=(2, 0))
        c.show()
    def compile(self):
        self.compile_electrodes()
        
    @abstractmethod
    def compile_electrodes(self):
        pass
