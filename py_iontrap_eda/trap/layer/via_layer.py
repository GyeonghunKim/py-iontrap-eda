from .base_layer import BaseLayer, LayerName
from typing import Dict, List
from ..electrode import Electrode
import shapely
from itertools import product
import numpy as np
class ViaLayer(BaseLayer):
    def __init__(self, name: LayerName, long_name: str, description: str, via_square_width: float, via_square_height: float, via_square_spacing: float, via_outside_margin: float):
        super().__init__(name, long_name, description)
        self.via_areas: List[shapely.Polygon] = []
        self.via_squares: List[shapely.Polygon] = []
        self.via_square_width = via_square_width
        self.via_square_height = via_square_height
        self.via_square_spacing = via_square_spacing
        self.via_outside_margin = via_outside_margin
        

        
    def add_via_area(self, via_area: shapely.Polygon):
        self.via_areas.append(via_area)
        
    def add_via_areas(self, via_areas: List[shapely.Polygon]):
        self.via_areas = self.via_areas + via_areas
        
    def __str__(self):
        return f"ViaLayer(name={self.name})"
    
    def __repr__(self):
        return f"ViaLayer(name={self.name}, via_areas={self.via_areas})"
    
    def arange_via_squares(self):
        for via_area in self.via_areas:
            min_x, min_y, max_x, max_y = via_area.bounds
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
            effective_width = max_x - min_x - 2 * self.via_outside_margin
            effective_height = max_y - min_y - 2 * self.via_outside_margin
            n_via_squares_width = int(effective_width / (self.via_square_width + self.via_square_spacing))
            n_via_squares_height = int(effective_height / (self.via_square_height + self.via_square_spacing))
            if n_via_squares_width % 2 == 0:
                x_centers = center_x + np.linspace(-(n_via_squares_width - 1)/2, (n_via_squares_width - 1)/2, n_via_squares_width) * (self.via_square_width + self.via_square_spacing)
            else:
                x_centers = center_x + np.linspace(-(n_via_squares_width - 1)/2, (n_via_squares_width - 1)/2, n_via_squares_width) * (self.via_square_width + self.via_square_spacing)
            
            if n_via_squares_height % 2 == 0:
                y_centers = center_y + np.linspace(-(n_via_squares_height - 1)/2, (n_via_squares_height - 1)/2, n_via_squares_height) * (self.via_square_height + self.via_square_spacing)
            else:
                y_centers = center_y + np.linspace(-(n_via_squares_height - 1)/2, (n_via_squares_height - 1)/2, n_via_squares_height) * (self.via_square_height + self.via_square_spacing)
            
            for x_center, y_center in product(x_centers, y_centers):
                via_square = shapely.Polygon([(x_center - self.via_square_width/2, y_center - self.via_square_height/2), (x_center + self.via_square_width/2, y_center - self.via_square_height/2), (x_center + self.via_square_width/2, y_center + self.via_square_height/2), (x_center - self.via_square_width/2, y_center + self.via_square_height/2)])
                self.via_squares.append(via_square)

    def arange_via_squares_in_rf_geometry(self, rf_geometry: shapely.Polygon):
        min_x, min_y, max_x, max_y = rf_geometry.bounds
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        effective_width = max_x - min_x - 2 * self.via_outside_margin
        effective_height = max_y - min_y - 2 * self.via_outside_margin
        n_via_squares_width = int(effective_width / (self.via_square_width + self.via_square_spacing))
        n_via_squares_height = int(effective_height / (self.via_square_height + self.via_square_spacing))
        bbox_via_squares = []
        if n_via_squares_width % 2 == 0:
            x_centers = center_x + np.linspace(-(n_via_squares_width - 1)/2, (n_via_squares_width - 1)/2, n_via_squares_width) * (self.via_square_width + self.via_square_spacing)
        else:
            x_centers = center_x + np.linspace(-(n_via_squares_width - 1)/2, (n_via_squares_width - 1)/2, n_via_squares_width) * (self.via_square_width + self.via_square_spacing)
        if n_via_squares_height % 2 == 0:
            y_centers = center_y + np.linspace(-(n_via_squares_height - 1)/2, (n_via_squares_height - 1)/2, n_via_squares_height) * (self.via_square_height + self.via_square_spacing)
        else:
            y_centers = center_y + np.linspace(-(n_via_squares_height - 1)/2, (n_via_squares_height - 1)/2, n_via_squares_height) * (self.via_square_height + self.via_square_spacing)
        for x_center, y_center in product(x_centers, y_centers):
            via_square = shapely.Polygon([(x_center - self.via_square_width/2, y_center - self.via_square_height/2), (x_center + self.via_square_width/2, y_center - self.via_square_height/2), (x_center + self.via_square_width/2, y_center + self.via_square_height/2), (x_center - self.via_square_width/2, y_center + self.via_square_height/2)])
            bbox_via_squares.append(via_square)
        
        rf_geometry_via_squares = []
        for square in bbox_via_squares:
            if rf_geometry.contains(square):
                rf_geometry_via_squares.append(square)
        return rf_geometry_via_squares