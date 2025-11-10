from ..trap.rail.base_rail import BaseRail
from ..utils import Units, Constants
import shapely
import numpy as np
from copy import deepcopy
from typing import Dict, List, Callable
from ..utils.ion_library import Ion
from functools import partial
from scipy.optimize import minimize
import numdifftools as nd

class GaplessROI:
    def __init__(self, x_min: float, x_max: float, y_min: float, y_max: float,):
        self.y_min = y_min
        self.y_max = y_max
        self.x_min = x_min
        self.x_max = x_max
        self.geometry = shapely.box(x_min, y_min, x_max, y_max)

class SquareElectrodeComponent:
    def __init__(self, x_min: float, x_max: float, y_min: float, y_max: float):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

class GaplessSolver:
    def __init__(self, rail: BaseRail, roi: GaplessROI, ion: Ion):
        self.rail = rail
        self.roi = roi
        self.electrodes_in_roi: Dict[str, List[SquareElectrodeComponent]] = {}
        self.get_intersection()
        self.potential_function_components: Dict[str, List[Callable[[float, float, float], float]]] = {}
        self.ion = ion
        self.get_potential_functions()
        self.get_unit_ponderomotive_potential()
        # self.ion_position = self.get_ion_position()
        # self.exit_position = self.get_exit_position()
        # self.get_unit_mathieu_Q()
        # self.get_unit_secular_frequency()
        # self.get_unit_trap_depth()
    def rf_potential_factory(self, rf_voltage: float):
        def rf_potential(x, y, z):
            eV = Constants.e * rf_voltage
            f_elec = self.potential_function('RF')
            return eV * f_elec(x, y, z)
        return rf_potential

    def potential_factory(self, rf_frequency: float, electrode_voltages: Dict[str, float]):

        def potential(x, y, z):
            terms = []

            for electrode_name, electrode_voltage in electrode_voltages.items():
                if electrode_name == 'RF':
                    rf_voltage = electrode_voltages['RF']
                    ponderomotive_potential = self.get_ponderomotive_potential(rf_voltage, rf_frequency)
                    terms.append(ponderomotive_potential)
                else:
                    eV = Constants.e * electrode_voltage
                    f_elec = self.potential_function(electrode_name)

                    terms.append(lambda x, y, z, eV=eV, f=f_elec: eV * f(x, y, z))

            return np.sum([t(x, y, z) for t in terms], axis=0)

        return potential
    
    # def get_unit_trap_depth(self):
    #     self.unit_trap_depth = (
    #         self.unit_ponderomotive_potential(0, 0, self.exit_position[0], self.exit_position[1]) - 
    #         self.unit_ponderomotive_potential(0, 0, self.ion_position[0], self.ion_position[1])
    #     )

    # def get_trap_depth(self, rf_voltage: float, rf_frequency: float):
    #     return rf_voltage**2 / rf_frequency**2 * self.unit_trap_depth

    # def scan_rf_voltage(self, rf_voltage_list: List[float], rf_frequency: float):
    #     trap_depth_list = []
    #     q_parameter_1_list = []
    #     q_parameter_2_list = []
    #     secular_frequency_1_list = []
    #     secular_frequency_2_list = []
    #     for rf_voltage in rf_voltage_list:
    #         trap_depth_list.append(self.get_trap_depth(rf_voltage, rf_frequency))
    #         q_1, q_2 = self.get_mathieu_Q(rf_voltage, rf_frequency)[:2]
    #         q_parameter_1_list.append(q_1)
    #         q_parameter_2_list.append(q_2)
    #         w1, w2 = self.get_secular_frequency(rf_voltage, rf_frequency)[:2]
    #         secular_frequency_1_list.append(w1)
    #         secular_frequency_2_list.append(w2)
    #     return np.array(trap_depth_list), np.array(q_parameter_1_list), np.array(q_parameter_2_list), np.array(secular_frequency_1_list), np.array(secular_frequency_2_list)
    
    # def get_axis_angle(self):
    #     return np.rad2deg(np.arccos(np.dot(self.unit_v1, [1, 0, 0]) / np.sqrt(np.dot(self.unit_v1, self.unit_v1))))
    
    # def get_unit_mathieu_Q(self):
    #     rf_potential_function = self.potential_function("RF")
    #     Q_matrix = 2 * self.ion.charge / (self.ion.ion_mass * (2 * np.pi * 1 * Units.Hz)**2) * nd.Hessian(lambda xy:rf_potential_function(xy[0], xy[1], 0), step=nd.step_generators.MaxStepGenerator(0.1*Units.um))(self.ion_position)
    #     eig_val, eig_vec = np.linalg.eig(Q_matrix)
    #     self.unit_Q = Q_matrix
    #     self.unit_q1 = eig_val[0]
    #     self.unit_q2 = eig_val[1]
    #     self.unit_v1 = eig_vec[:,0]
    #     self.unit_v2 = eig_vec[:,1]
    
    # def get_mathieu_Q(self, rf_voltage: float, rf_frequency: float):
    #     Q = rf_voltage / rf_frequency**2 * self.unit_Q
    #     q1 = rf_voltage / rf_frequency**2 * self.unit_q1
    #     q2 = rf_voltage / rf_frequency**2 * self.unit_q2
    #     return Q, q1, q2, self.unit_v1, self.unit_v2

    # def get_unit_secular_frequency(self):
    #     hessian = nd.Hessian(lambda xy:self.unit_ponderomotive_potential(xy[0], xy[1], 0), step=nd.step_generators.MaxStepGenerator(0.1*Units.um))(self.get_ion_position())
    #     eig_val, eig_vec = np.linalg.eig(hessian)
    #     sec_freq = np.sqrt(np.abs(eig_val[:2]) / self.ion.ion_mass)
    #     self.unit_secular_frequency_1 = sec_freq[0]
    #     self.unit_secular_frequency_2 = sec_freq[1]
    #     self.unit_secular_axis_1 = eig_vec[:,0]
    #     self.unit_secular_axis_2 = eig_vec[:,1]

    # def get_secular_frequency(self, rf_voltage: float, rf_frequency: float):
    #     return rf_voltage / rf_frequency * self.unit_secular_frequency_1, rf_voltage / rf_frequency * self.unit_secular_frequency_2, self.unit_secular_axis_1, self.unit_secular_axis_2
    
    # def get_exit_position(self, dx=100*Units.nm, dy=100*Units.nm, dz=100*Units.nm):
    #     ponderomotive_potential = self.get_ponderomotive_potential(rf_voltage=100*Units.V, rf_frequency=35*Units.MHz)
    #     def cost(xy):
    #         dpsi_dx = (ponderomotive_potential(xy[0] + dx, xy[1], 0) - ponderomotive_potential(xy[0] - dx, xy[1], 0)) / (2 * dx)
    #         dpsi_dy = (ponderomotive_potential(xy[0], xy[1] + dy, 0) - ponderomotive_potential(xy[0], xy[1] - dy, 0)) / (2 * dy)
    #         return (dpsi_dx**2 + dpsi_dy**2)
    #     res = minimize(cost, [2 * Units.um, 130 * Units.um], method="Nelder-Mead", options={"xatol": 1e-35, "disp": False})
    #     exit_position = res.x
    #     return np.array([exit_position[0], exit_position[1], 0])
    
    # def get_ion_position(self, dx=100*Units.nm, dy=100*Units.nm, dz=100*Units.nm):
    #     def cost(xy):
    #         dpsi_dx = (self.unit_ponderomotive_potential(xy[0] + dx, xy[1], 0) - self.unit_ponderomotive_potential(xy[0] - dx, xy[1], 0)) / (2 * dx)
    #         dpsi_dy = (self.unit_ponderomotive_potential(xy[0], xy[1] + dy, 0) - self.unit_ponderomotive_potential(xy[0], xy[1] - dy, 0)) / (2 * dy)
    #         return dpsi_dx**2 + dpsi_dy**2
    #     res = minimize(cost, [1*Units.um, 70 * Units.um], method="Nelder-Mead", options={"xatol": 1e-35, "disp": False})
    #     ion_position = res.x
    #     return np.array([ion_position[0], ion_position[1], 0])
    
    def get_unit_ponderomotive_potential(self, dx=100*Units.nm, dy=100*Units.nm, dz=100*Units.nm):
        rf_potential_function = self.potential_function("RF")

        def ponderomotive_potential(x, y, z):
            dPhi_dx = (rf_potential_function(x + dx, y, z) - rf_potential_function(x - dx, y, z)) / (2 * dx)
            dPhi_dy = (rf_potential_function(x, y + dy, z) - rf_potential_function(x, y - dy, z)) / (2 * dy)
            dPhi_dz = (rf_potential_function(x, y, z + dz) - rf_potential_function(x, y, z - dz)) / (2 * dz)
            return self.ion.charge**2 * 1*Units.V**2 / (4 * self.ion.ion_mass * (2 * np.pi * 1 * Units.Hz) ** 2) * (dPhi_dx**2 + dPhi_dy**2 + dPhi_dz**2)
        self.unit_ponderomotive_potential = ponderomotive_potential

    def get_ponderomotive_potential(self, rf_voltage: float, rf_frequency: float):
        ponderomotive_potential = lambda x, y, z: self.unit_ponderomotive_potential(x, y, z) * rf_voltage**2 / (rf_frequency**2)
        return ponderomotive_potential


    def get_intersection(self):
        for electrode in self.rail.electrodes:
            self.electrodes_in_roi[electrode.name] = []
            intersection = self.roi.geometry.intersection(electrode.geometry)
            if intersection.area > 1e-20:   
                if isinstance(intersection, shapely.MultiPolygon):
                    for poly in intersection.geoms:
                        x_min, y_min, x_max, y_max = poly.bounds
                        if np.isclose(poly.area, (x_max - x_min) * (y_max - y_min)):
                            self.electrodes_in_roi[electrode.name].append(
                                SquareElectrodeComponent(x_min, x_max, y_min, y_max)
                            )
                else:
                    x_min, y_min, x_max, y_max = intersection.bounds
                    if np.isclose(intersection.area, (x_max - x_min) * (y_max - y_min)):
                        self.electrodes_in_roi[electrode.name].append(
                            SquareElectrodeComponent(x_min, x_max, y_min, y_max)
                        )
            else:   
                self.electrodes_in_roi[electrode.name].append(None)

    def get_potential_functions(self):
        for electrode_name, square_electrode_components in self.electrodes_in_roi.items():
            self.potential_function_components[electrode_name] = []
            for square_electrode_component in square_electrode_components:
                self.potential_function_components[electrode_name].append(self.get_component_solution(square_electrode_component))
            
    def potential_function(self, electrode_name: str):
        potential_functions = deepcopy(self.potential_function_components[electrode_name])
        def potential_function(x, y, z):
            return np.sum([potential_function_component(x, y, z) for potential_function_component in potential_functions], axis=0)
        return potential_function

    def get_component_solution(self, square_electrode_component: SquareElectrodeComponent):

        x1 = square_electrode_component.x_min
        x2 = square_electrode_component.x_max
        y1 = square_electrode_component.y_min
        y2 = square_electrode_component.y_max
        def solution(x, y, z):
            term_1 = np.arctan2((y2 - y) * (x2 - x), z * np.sqrt(z**2 + (y2 - y) ** 2 + (x2 - x) ** 2))
            term_2 = np.arctan2((y1 - y) * (x2 - x), z * np.sqrt(z**2 + (y1 - y) ** 2 + (x2 - x) ** 2))
            term_3 = np.arctan2((y2 - y) * (x1 - x), z * np.sqrt(z**2 + (y2 - y) ** 2 + (x1 - x) ** 2))
            term_4 = np.arctan2((y1 - y) * (x1 - x), z * np.sqrt(z**2 + (y1 - y) ** 2 + (x1 - x) ** 2))
            return (term_1 - term_2 - term_3 + term_4) / (2 * np.pi)
        return solution
