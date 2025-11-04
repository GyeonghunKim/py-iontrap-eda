from dataclasses import dataclass
from typing import List

import numpy as np
import shapely

from .base_rail import BaseRail, BaseRailParameters
from ..electrode import Electrode
from .zone import Zone, Spacing
from ...utils import Units


@dataclass
class TwoInnerDCRailParameters(BaseRailParameters):
    n_layers: int
    inner_dc_width: float
    rf_width: float
    rf_merge_length: float
    rf_padding_left_width: float
    rf_padding_right_width: float


class TwoInnerDCRail(BaseRail):
    def __init__(self, name: str, parameters: TwoInnerDCRailParameters):
        self.parameters: TwoInnerDCRailParameters
        super().__init__(name, parameters)

    def compile_electrodes(self):
        self.electrodes = []

        # RF
        self.electrodes.append(
            Electrode(
                "RF",
                shapely.Polygon(
                    (
                        (
                            self.dc_z_points[0] - self.parameters.rf_padding_left_width - self.parameters.rf_merge_length,
                            self.parameters.inner_dc_width + self.parameters.rf_width,
                        ),
                        (
                            self.dc_z_points[-1] + self.parameters.rf_padding_right_width   ,
                            self.parameters.inner_dc_width + self.parameters.rf_width,
                        ),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, self.parameters.inner_dc_width),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, self.parameters.inner_dc_width),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, -self.parameters.inner_dc_width),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, -self.parameters.inner_dc_width),
                        (
                            self.dc_z_points[-1] + self.parameters.rf_padding_right_width,
                            -self.parameters.inner_dc_width - self.parameters.rf_width,
                        ),
                        (
                            self.dc_z_points[0] - self.parameters.rf_padding_left_width - self.parameters.rf_merge_length,
                            -self.parameters.inner_dc_width - self.parameters.rf_width,
                        ),
                        (
                            self.dc_z_points[0] - self.parameters.rf_padding_left_width - self.parameters.rf_merge_length,
                            self.parameters.inner_dc_width + self.parameters.rf_width,
                        ),
                    )
                ),
            )
        )

        # INNER DC
        self.electrodes.append(
            Electrode(
                "E0",
                shapely.Polygon(
                    (
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, self.parameters.inner_dc_width),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, self.parameters.inner_dc_width),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, self.parameters.inner_dc_width + self.parameters.rf_width),
                        (self.dc_z_points[-1] + self.parameters.z_padding + self.parameters.rf_padding_right_width, self.parameters.inner_dc_width + self.parameters.rf_width),
                        (self.dc_z_points[-1] + self.parameters.z_padding + self.parameters.rf_padding_right_width, 0),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, 0),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, self.parameters.inner_dc_width),
                    )
                ),
            )
        )
        self.electrodes.append(
            Electrode(
                "E1",
                shapely.Polygon(
                    (
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, -self.parameters.inner_dc_width),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, -self.parameters.inner_dc_width),
                        (self.dc_z_points[-1] + self.parameters.rf_padding_right_width, -(self.parameters.inner_dc_width + self.parameters.rf_width)),
                        (self.dc_z_points[-1] + self.parameters.z_padding + self.parameters.rf_padding_right_width, -(self.parameters.inner_dc_width + self.parameters.rf_width)),
                        (self.dc_z_points[-1] + self.parameters.z_padding + self.parameters.rf_padding_right_width, 0),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, 0),
                        (self.dc_z_points[0] - self.parameters.rf_padding_left_width, -self.parameters.inner_dc_width),
                    )
                ),
            )
        )

        # OUTER DC
        electrode_count = 2
        for z_left, z_right in zip(self.dc_z_points[:-1], self.dc_z_points[1:]):
            for sign in [1, -1]:
                self.electrodes.append(
                    Electrode(
                        f"E{electrode_count}",
                        shapely.Polygon(
                            (
                                (z_left, sign * (self.parameters.inner_dc_width + self.parameters.rf_width + self.parameters.dc_length)),
                                (z_right, sign * (self.parameters.inner_dc_width + self.parameters.rf_width + self.parameters.dc_length)),
                                (z_right, sign * (self.parameters.inner_dc_width + self.parameters.rf_width)),
                                (z_left, sign * (self.parameters.inner_dc_width + self.parameters.rf_width)),
                                (z_left, sign * (self.parameters.inner_dc_width + self.parameters.rf_width + self.parameters.dc_length)),
                            )
                        ),
                    )
                )
                electrode_count += 1
                
        self.generate_ports()
        if self.parameters.via_pad_width is not None and self.parameters.via_pad_gap is not None:   
            self.generate_via_pads()
