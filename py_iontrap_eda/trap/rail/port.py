import shapely

class Port:
    def __init__(self, name: str, geometry: shapely.LineString, orientation: int): # orientation 0 for right, 180 for left, 90 for up, -90 for down
        self.name = name
        self.geometry = geometry
        self.orientation = orientation
        x_min, y_min, x_max, y_max = self.geometry.bounds
        self.center = (x_min + x_max) / 2, (y_min + y_max) / 2
        self.width = abs(y_max - y_min) if orientation == 0 or orientation == 180 else abs(x_max - x_min)
        
    def __str__(self):
        return f"{self.name}: {self.geometry}"

    def __repr__(self):
        return f"Port(name={self.name}, geometry={self.geometry}, orientation={self.orientation})"

    def __eq__(self, other):
        return self.name == other.name and self.geometry == other.geometry and self.orientation == other.orientation

    def __hash__(self):
        return hash((self.name, self.geometry, self.orientation))