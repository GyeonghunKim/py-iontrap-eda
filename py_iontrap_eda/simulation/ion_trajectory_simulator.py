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

    def get_dc_electric_force(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, eps: float = 0.1*Units.um):
        x_list = [x - eps, x + eps, x, x, x, x]
        y_list = [y, y, y - eps, y + eps, y, y]
        z_list = [z, z, z, z, z - eps, z + eps]
        potential_list = [self.dc_potential_function(x_i, y_i, z_i) for x_i, y_i, z_i in zip(x_list, y_list, z_list)]
        force_x = (potential_list[0] - potential_list[1]) / (2 * eps)
        force_y = (potential_list[2] - potential_list[3]) / (2 * eps)
        force_z = (potential_list[4] - potential_list[5]) / (2 * eps)
        return np.array([force_x, force_y, force_z])

    def get_rf_electric_force(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, eps: float = 0.1*Units.um):
        x_list = [x - eps, x + eps, x, x, x, x]
        y_list = [y, y, y - eps, y + eps, y, y]
        z_list = [z, z, z, z, z - eps, z + eps]
        potential_list = [self.rf_potential_function(x_i, y_i, z_i) for x_i, y_i, z_i in zip(x_list, y_list, z_list)]
        force_x = (potential_list[0] - potential_list[1]) / (2 * eps)
        force_y = (potential_list[2] - potential_list[3]) / (2 * eps)
        force_z = (potential_list[4] - potential_list[5]) / (2 * eps)
        return np.array([force_x, force_y, force_z])

    def solve_exact_custom_rk4(self, x0_list: np.ndarray, y0_list: np.ndarray, z0_list: np.ndarray, vx0_list: np.ndarray, vy0_list: np.ndarray, vz0_list: np.ndarray, t_list: np.ndarray):
        phase_space_initial_conditions = np.stack([x0_list, y0_list, z0_list, vx0_list, vy0_list, vz0_list], axis=1)
        def rk4_step(fun, t, y, dt):
            k1 = fun(t, y)
            k2 = fun(t + dt/2, y + dt/2 * k1)
            k3 = fun(t + dt/2, y + dt/2 * k2)
            k4 = fun(t + dt,   y + dt * k3)
            return y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

        def f(t, phase_space_variable):
            x, y, z= phase_space_variable[:, :3].T
            vx, vy, vz = phase_space_variable[:, 3:].T
            dxdt = vx
            dydt = vy
            dzdt = vz
            dvdt = -(self.get_dc_electric_force(x, y, z) + self.get_rf_electric_force(x, y, z) * np.cos(2 * np.pi * self.rf_frequency * t)) / (self.ion.ion_mass)
            return np.array([dxdt, dydt, dzdt, *dvdt]).T
        res = [phase_space_initial_conditions]
        for i, t in enumerate(t_list[1:]):
            res.append(rk4_step(f, t, res[-1], t_list[1] - t_list[0]))

        return np.array(res)

    def solve_exact(self, initial_position: np.ndarray, initial_velocity: np.ndarray, t_list: np.ndarray):
        def f(t, phase_space_variable):
            x, y, z= phase_space_variable[:3]
            r = np.array([x, y, z])
            vx, vy, vz = phase_space_variable[3:]
            dxdt = vx
            dydt = vy
            dzdt = vz
            dvdt = -(self.get_dc_electric_force(x, y, z) + self.get_rf_electric_force(x, y, z) * np.cos(2 * np.pi * self.rf_frequency * t)) / (self.ion.ion_mass)
            return np.array([dxdt, dydt, dzdt, *dvdt.flatten()])
        res = solve_ivp(f, [t_list[0], t_list[-1]], np.hstack([initial_position, initial_velocity]), method='RK45', t_eval=t_list)
        return res

    def get_stiffness_matrix(self):
        prefactor = 4 / (self.ion.ion_mass * (2 * np.pi * self.rf_frequency)**2)
        hessian = nd.Hessian(lambda xyz:self.dc_potential_function(xyz[0], xyz[1], xyz[2]), step=nd.step_generators.MaxStepGenerator(2*Units.um))
        return lambda x, y, z: prefactor * hessian(np.array([x, y, z]))

    def get_excitation_matrix(self):
        prefactor = 2 / (self.ion.ion_mass * (2 * np.pi * self.rf_frequency)**2)
        hessian = nd.Hessian(lambda xyz:self.rf_potential_function(xyz[0], xyz[1], xyz[2]), step=nd.step_generators.MaxStepGenerator(2*Units.um))

        return lambda x, y, z: prefactor *  hessian(np.array([x, y, z]))
    