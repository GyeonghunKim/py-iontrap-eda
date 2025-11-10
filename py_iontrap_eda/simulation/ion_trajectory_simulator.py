from typing import Callable, Tuple
from ..utils import Units, Constants
from ..utils.ion_library import Ion
import numpy as np
import numdifftools as nd
from scipy.integrate import solve_ivp
class IonTrajectorySimulator:
    def __init__(self, 
    rf_potential_function: Callable[[float, float, float], float], 
    dc_potential_function: Callable[[float, float, float], float],
    ion: Ion,
    rf_frequency: float,
    ):
        self.rf_potential_function = rf_potential_function
        self.dc_potential_function = dc_potential_function
        self.ion = ion
        self.rf_frequency = rf_frequency
        
        self.stiffness_matrix_function = self.get_stiffness_matrix()
        self.excitation_matrix_function = self.get_excitation_matrix()

    def get_stiffness_matrix(self):
        prefactor = 4 * self.ion.charge / (self.ion.ion_mass * (2 * np.pi * self.rf_frequency)**2)
        hessian = nd.Hessian(lambda xyz:self.dc_potential_function(xyz[0], xyz[1], xyz[2]), step=nd.step_generators.MaxStepGenerator(1*Units.um))
        return lambda x, y, z: prefactor * hessian((x, y, z))

    def get_excitation_matrix(self):
        prefactor = 2 * self.ion.charge / (self.ion.ion_mass * (2 * np.pi * self.rf_frequency)**2)
        hessian = nd.Hessian(lambda xyz:self.rf_potential_function(xyz[0], xyz[1], xyz[2]), step=nd.step_generators.MaxStepGenerator(1*Units.um))
        return lambda x, y, z: prefactor *  hessian((x, y, z))

