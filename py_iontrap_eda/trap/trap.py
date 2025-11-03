from typing import List, Dict, Optional, Tuple, Callable
from .layer.base_layer import BaseLayer, LayerName
from .thickness_map import ThicknessMap
from .rail.base_rail import BaseRail
import shapely

import gdsfactory as gf
import numpy as np

from ..utils import Units
from .bond_pad import BondPad
from .route_method import RouteMethod
from .electrode import Electrode
class Trap:
    def __init__(self, 
                 name: str, 
                 description: str, 
                 width: float,
                 height: float,
                 dicing_margin: float,
                 dicing_channel_layer: LayerName,
                 layers: Dict[LayerName, BaseLayer], 
                 thickness_map: ThicknessMap):
        self.name = name
        self.description = description
        self.width = width
        self.height = height
        self.dicing_margin = dicing_margin
        self.layers = layers
        self.thickness_map = thickness_map  
        self.dicing_channel_layer = dicing_channel_layer
        self.create_dicing_channel()
        self.bond_pads: Dict[str, BondPad] = {}
        
    def create_dicing_channel(self):
        self.die_geometry = shapely.Polygon([
            (-self.width/2, -self.height/2),
            (self.width/2, -self.height/2),
            (self.width/2, self.height/2),
            (-self.width/2, self.height/2),
        ])
        self.dicing_channel_geometry = shapely.Polygon(
            (
                (-self.width/2 - self.dicing_margin, -self.height/2 - self.dicing_margin),
                (self.width/2 + self.dicing_margin, -self.height/2 - self.dicing_margin),
                (self.width/2 + self.dicing_margin, self.height/2 + self.dicing_margin),
                (-self.width/2 - self.dicing_margin, self.height/2 + self.dicing_margin),
            )
        )
        self.dicing_channel_geometry = self.dicing_channel_geometry.difference(self.die_geometry)
        
    def add_rail(self, rail: BaseRail):
        self.rail = rail
        top_layer_name = self.thickness_map.layer_names[0]
        top_layer = self.layers[top_layer_name]
        top_layer.add_total_negatives(list(rail.total_negatives.geoms))
        for electrode in rail.electrodes:
            top_layer.add_electrode(
                electrode
            )
    def add_bond_pad(self, bond_pad: BondPad):
        self.bond_pads[bond_pad.name] = bond_pad
        electrode = self.rail.get_electrode(bond_pad.name)
        if electrode is not None and electrode.route_method == RouteMethod.VIA:
            electrode.via_pad_layer_down_to = bond_pad.via_pad_layer_down_to
            
            
    def add_bond_pads(self, 
                      names: List[str], 
                      routing_methods: List[RouteMethod], 
                      location: str, 
                      width: float, 
                      height: float, 
                      offset_x: float, 
                      offset_y: float, 
                      gap: float,
                      via_pad_centers: Optional[List[Tuple[float, float]]]=None,
                      via_pad_widths: Optional[List[float]]=None,
                      via_pad_heights: Optional[List[float]]=None,
                      via_pad_layer_down_tos: Optional[List[LayerName]]=None,
                      ):
        if via_pad_centers is None or via_pad_widths is None or via_pad_heights is None or via_pad_layers is None:
            via_pad_centers = [None] * len(names)
            via_pad_widths = [None] * len(names)
            via_pad_heights = [None] * len(names)
            via_pad_layers = [None] * len(names)
        dr = np.array([0, 0])
        
        if location == "left-top":
            x = -self.width/2 + offset_x + width/2
            y = 0 + offset_y + height/2
            dr = np.array([0, width + gap])
            orientation = 0
            
        elif location == "left-bottom":
            x = -self.width/2 + offset_x + width/2
            y = 0 - offset_y - height/2
            dr = np.array([0, - (width + gap)])
            orientation = 0
            
        else:
            raise ValueError(f"Invalid location: {location}")
        
        for i, (name, routing_method, via_pad_center, via_pad_width, via_pad_height, via_pad_layer_down_to) in enumerate(
            zip(names, routing_methods, via_pad_centers, via_pad_widths, via_pad_heights, via_pad_layer_down_tos)
            ):
            if via_pad_center is None:
                via_pad_center = (x + i * dr[0], y + i * dr[1])
            if via_pad_width is None:
                via_pad_width = width
            if via_pad_height is None:
                via_pad_height = height
            
            self.add_bond_pad(BondPad(name, (x + i * dr[0], y + i * dr[1]), width, height, orientation, routing_method, via_pad_center, via_pad_width, via_pad_height, via_pad_layer_down_to))
            if routing_method == RouteMethod.VIA:
                electrode = self.rail.get_electrode(name)
                if electrode is not None:
                    electrode.via_pad_layer_down_to = via_pad_layer_down_to
    def compile_via_pads(self):
        for electrode in self.rail.electrodes:
            if electrode.route_method == RouteMethod.VIA:
                self.layers[electrode.via_pad_layer_down_to].add_electrode(Electrode(electrode.name, electrode.via_pad, RouteMethod.VIA))
                
    def compile_dicing_channel(self):
        def dicing_channel_factory():
            dicing_channel = gf.Component()
            dicing_channel_outer = gf.Component()
            dicing_channel_inner = gf.Component()
            dicing_channel_outer.add_polygon(np.array(self.dicing_channel_geometry.exterior.xy).T / Units.um, layer=self.dicing_channel_layer)
            dicing_channel_inner.add_polygon(np.array(self.die_geometry.exterior.xy).T / Units.um, layer=self.dicing_channel_layer)
            dicing_channel = gf.boolean(dicing_channel_outer, dicing_channel_inner, 'A-B', layer=self.dicing_channel_layer)
            return dicing_channel
        
        self.dicing_channel_factory = dicing_channel_factory
    
    def add_auto_port_routing(self, electrode_name: str, bond_pad_name: str):
        electrode = self.rail.get_electrode(electrode_name)
        electrode_component = gf.Component()
        electrode_component.add_polygon(np.array(electrode.geometry.exterior.xy).T / Units.um, layer = (780, 727))
        electrode_component.add_port(name=electrode_name, 
                    center=np.array(electrode.port.center) / Units.um, 
                    width=electrode.port.width / Units.um, 
                    orientation=electrode.port.orientation,
                    layer=(780, 727)
                )
        pad_component = gf.Component()
        pad_component.add_polygon(np.array(self.bond_pads[bond_pad_name].geometry.exterior.xy).T / Units.um, layer = (780, 727))
        pad_component.add_port(name=bond_pad_name, 
                center=np.array(self.bond_pads[bond_pad_name].port.center) / Units.um, 
                width=self.bond_pads[bond_pad_name].port.width / Units.um, 
                orientation=self.bond_pads[bond_pad_name].port.orientation,
                layer=(780, 727)
                )
                
        c = gf.Component()
        c << pad_component
        c << electrode_component

        route = gf.routing.route_single_electrical(
            c,
            pad_component.ports[bond_pad_name],
            electrode_component.ports[electrode_name],
            # cross_section=gf.cross_section.strip,
            width=pad_component.ports[bond_pad_name].width,
            layer=(780, 727)
        )
        c.show()
        shapely_polygons = []
        for polygons in list(c.get_polygons_points().values()):
            for polygon in polygons:
                shapely_polygons.append(shapely.Polygon(polygon * Units.um))
        total_polygon = shapely.union_all(shapely_polygons)
        self.layers[self.thickness_map.layer_names[0]].add_electrode(Electrode(electrode_name, total_polygon, RouteMethod.ROUTED))
        self.bond_pads.pop(bond_pad_name)
        
    def compile_bond_pads(self):
        top_layer_name = self.thickness_map.layer_names[0]
        def bond_pad_factory(bond_pad: BondPad):
            c = gf.Component()
            c.add_polygon(np.array(bond_pad.geometry.exterior.xy).T / Units.um, layer=top_layer_name)
            c.add_label(bond_pad.name, position=c.dcenter, layer=top_layer_name)
            if bond_pad.route_method == RouteMethod.PORT:
                c.add_port(name=bond_pad.name, center=np.array(bond_pad.port.center) / Units.um, width=bond_pad.port.width / Units.um, orientation=bond_pad.port.orientation, layer=top_layer_name)
            return c
        
        def bond_pads_factory():
            c = gf.Component()
            for bond_pad in self.bond_pads.values():
                c_temp = bond_pad_factory(bond_pad)
                c << c_temp
            return c
        
        self.bond_pads_factory = bond_pads_factory
    def compile_rf_ground_gap_tapering(self, tapering_length_from_left_edge: float, max_gap: float, min_gap: float, max_gap_length: float):
        rf_pad_height =(self.rail.parameters.inner_dc_width + self.rail.parameters.rf_width)

        upper_polygon = shapely.Polygon(
            [(-self.width/2, rf_pad_height + max_gap - min_gap / 2), 
             (-self.width/2 + max_gap_length, rf_pad_height + max_gap - min_gap / 2), 
             (-self.width/2 + tapering_length_from_left_edge, rf_pad_height + min_gap / 2), 
             (-self.width/2 + tapering_length_from_left_edge, rf_pad_height - min_gap / 2), 
             (-self.width/2, rf_pad_height - min_gap / 2)]
        )
        lower_polygon = shapely.Polygon(
            [(-self.width/2, -(rf_pad_height + max_gap - min_gap / 2)), 
             (-self.width/2 + max_gap_length, -(rf_pad_height + max_gap - min_gap / 2)), 
             (-self.width/2 + tapering_length_from_left_edge, -(rf_pad_height + min_gap / 2)), 
             (-self.width/2 + tapering_length_from_left_edge, -(rf_pad_height - min_gap / 2)), 
             (-self.width/2, -(rf_pad_height - min_gap / 2))]
        )
        rf_ground_gap_tapering = shapely.union_all([upper_polygon, lower_polygon])
        self.top_layer_grounds = shapely.difference(self.top_layer_grounds, rf_ground_gap_tapering)
        top_layer_name = self.thickness_map.layer_names[0]
        top_layer = self.layers[top_layer_name]
        top_layer.electrodes["RF"].geometry = shapely.difference(top_layer.electrodes["RF"].geometry, rf_ground_gap_tapering)
        
    def compile_top_layer_grounds(self, 
                                  left_keep: Optional[float]=None, 
                                  right_keep: Optional[float]=None, 
                                  top_keep: Optional[float]=None, 
                                  bottom_keep: Optional[float]=None,
                                  keep_electrodes: Optional[List[str]]=None,
                                  custom_keep_polygons: Optional[List[shapely.Polygon]]=None,
                                  ):
        if custom_keep_polygons is None:
            custom_keep_polygons = []
        top_layer_name = self.thickness_map.layer_names[0]
        top_layer = self.layers[top_layer_name]
        die_geometry = self.die_geometry

            
        all_polygons = []
        for electrode in top_layer.electrodes.values():
            all_polygons.append(electrode.geometry)
        for bond_pad in self.bond_pads.values():
            all_polygons.append(bond_pad.geometry)

        total_polygon = shapely.union_all(all_polygons)
        top_layer_grounds = shapely.difference(die_geometry, total_polygon)
        if left_keep is not None:
            top_layer_grounds = top_layer_grounds.difference(
                shapely.Polygon(
                    [(-self.width/2, -self.height/2), 
                     (-self.width/2, self.height/2), 
                     (-self.width/2 + left_keep, self.height/2), 
                     (-self.width/2 + left_keep, -self.height/2)]
                )
            )
        if right_keep is not None:
            top_layer_grounds = top_layer_grounds.difference(
                shapely.Polygon(
                    [(self.width/2, -self.height/2), 
                     (self.width/2, self.height/2), 
                     (self.width/2 - right_keep, self.height/2), 
                     (self.width/2 - right_keep, -self.height/2)]
                )
            )
        if top_keep is not None:
            top_layer_grounds = top_layer_grounds.difference(
                shapely.Polygon(
                    [(-self.width/2, self.height/2), 
                     (self.width/2, self.height/2), 
                     (self.width/2, self.height/2 - top_keep), 
                     (-self.width/2, self.height/2 - top_keep)]
                )
            )
        if bottom_keep is not None:
            top_layer_grounds = top_layer_grounds.difference(
                shapely.Polygon(
                    [(-self.width/2, -self.height/2), 
                     (self.width/2, -self.height/2), 
                     (self.width/2, -self.height/2 + bottom_keep), 
                     (-self.width/2, -self.height/2 + bottom_keep)]
                )
            )
            
        for polygon in custom_keep_polygons:
            top_layer_grounds = top_layer_grounds.difference(polygon)
            
        self.top_layer_grounds = top_layer_grounds
            
            
    def compile_top_layer_gaps(self, gap_function: Callable[[shapely.Polygon], shapely.Polygon], keep_electrodes: Optional[List[str]]=None, gap_join_style: str='mitre'):
        top_layer_name = self.thickness_map.layer_names[0]
        top_layer = self.layers[top_layer_name]
        all_polygons = []
        for electrode in top_layer.electrodes.values():
            all_polygons.append(electrode.geometry)
        for bond_pad in self.bond_pads.values():
            all_polygons.append(bond_pad.geometry)
        all_polygons.append(self.top_layer_grounds)
        
        boundaries = []
        for i, one in enumerate(all_polygons[:-1]):
            for two in all_polygons[i+1:]:
                inter = shapely.intersection(one, two, grid_size=1e-10)
                if inter.length > 1e-16:
                    if isinstance(inter, shapely.MultiLineString):
                        boundaries += list(inter.geoms)
                    else:
                        boundaries.append(inter)
        # return boundaries
        if keep_electrodes is None:
            keep_electrodes = []
        
        gaps = shapely.buffer(
            boundaries,
            [gap_function(boundary) for boundary in boundaries],
            join_style=gap_join_style
        )
        
        total_boundary = shapely.union_all(gaps)
        # return total_boundary
        for electrode_name, electrode in top_layer.electrodes.items():
            if electrode_name in keep_electrodes:
                continue
            else:
                top_layer.electrodes[electrode_name].geometry = shapely.difference(electrode.geometry, total_boundary)
        self.top_layer_grounds = shapely.difference(self.top_layer_grounds, total_boundary)
    
    def compile_top_layer_factory(self):
        top_layer_name = self.thickness_map.layer_names[0]
        top_layer = self.layers[top_layer_name]
        def top_layer_factory():
            c = gf.Component()
            # print(self.top_layer_grounds)
            if isinstance(self.top_layer_grounds, shapely.MultiPolygon):
                for geometry in self.top_layer_grounds.geoms:
                    exterior = gf.Component()
                    exterior.add_polygon(np.array(geometry.exterior.xy).T / Units.um, layer=top_layer_name)
                    interior = gf.Component()
                    for one_interior in geometry.interiors:
                        interior_component = gf.Component()
                        interior_component.add_polygon(np.array(one_interior.xy).T / Units.um, layer=top_layer_name)
                        interior = gf.boolean(interior, interior_component, 'or', layer=top_layer_name)
                    c << gf.boolean(exterior, interior, 'A-B', layer=top_layer_name)
            else:
                exterior = gf.Component()
                exterior.add_polygon(np.array(self.top_layer_grounds.exterior.xy).T / Units.um, layer=top_layer_name)
                interior = gf.Component()
                for one_interior in self.top_layer_grounds.interiors:
                    interior_component = gf.Component()
                    interior_component.add_polygon(np.array(one_interior.xy).T / Units.um, layer=top_layer_name)
                    interior = gf.boolean(interior, interior_component, 'or', layer=top_layer_name)
                c << gf.boolean(exterior, interior, 'A-B', layer=top_layer_name)
            for electrode in top_layer.electrodes.values():
                c.add_polygon(np.array(electrode.geometry.exterior.xy).T / Units.um, layer=top_layer_name)
            for bond_pad in self.bond_pads.values():
                c.add_polygon(np.array(bond_pad.geometry.exterior.xy).T / Units.um, layer=top_layer_name)
            return c
        self.top_layer_factory = top_layer_factory
        
        
    def compile(self):
        self.compile_dicing_channel()
        self.compile_top_layer_factory()
        
        def trap_factory():
            c = gf.Component()
            c << self.dicing_channel_factory()
            c << self.top_layer_factory()
            return c
        
        self.trap_factory = trap_factory
        
    def write_gds(self, file_name: str):        
        c = self.trap_factory()
        c.show()
        # c.write_gds(file_name)