import numpy as np
import gdstk

from .units import Units


def text_gds(text: str, height: float, center_x: float, center_y: float):
    text = gdstk.text(text, height / Units.um, (0, 0))
    min_x = 1e10
    min_y = 1e10
    max_x = -1e10
    max_y = -1e10
    for polygon in text:
        bbox = polygon.bounding_box()
        min_x = min(min_x, bbox[0][0])
        min_y = min(min_y, bbox[0][1])
        max_x = max(max_x, bbox[1][0])
        max_y = max(max_y, bbox[1][1])
        
    length = max_x - min_x
    height = max_y - min_y
    
    new_text = []
    for polygon in text:
        new_text.append(polygon.translate(-min_x - length / 2 + center_x / Units.um, -min_y - height / 2 + center_y / Units.um))
    return new_text
def duke_logo_gds(width: float, center_x: float, center_y: float):
    scale = width / Units.um / 270.001
    x_offset = -270.001/2
    y_offset = -240.223/2
    def transform_x(x):
        return (x + x_offset) * scale + center_x / Units.um
    def transform_y(y):
        return (y + y_offset) * scale + center_y / Units.um
    left_body = gdstk.Polygon(np.array([
        [transform_x(145), transform_y(320.108)],
        [transform_x(287.79), transform_y(320.108)],
        [transform_x(254.018), transform_y(264.446)],
        [transform_x(254.018), transform_y(135.719)],
        [transform_x(287.79), transform_y(79.918)],
        [transform_x(145.), transform_y(79.918)],
        [transform_x(178.759), transform_y(135.719)],
        [transform_x(178.759), transform_y(264.446)]
    ]))
    right_body = gdstk.Curve((transform_x(328.945), transform_y(264.446)), tolerance=1e-3)
    right_body.commands(
        "L", transform_x(295.221), transform_y(320.115), 
        "L", transform_x(365.868), transform_y(320.115), 
        "C", transform_x(365.868), transform_y(320.115), transform_x(400.951), transform_y(304.754), transform_x(415.001), transform_y(268.731),
        "L", transform_x(415.001), transform_y(130.217),
        "C", transform_x(415.001), transform_y(130.217), transform_x(400.956), transform_y(94.739), transform_x(365.868), transform_y(79.885),
        "L", transform_x(295.221), transform_y(79.885),
        "L", transform_x(328.945), transform_y(135.719),
    )
    
    right_body = gdstk.Polygon(right_body.points())
    return left_body, right_body

