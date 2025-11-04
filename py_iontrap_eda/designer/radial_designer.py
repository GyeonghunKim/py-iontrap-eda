
from typing import Optional, Tuple, List
from enum import Enum, auto
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from .design import RadialDesign
from ..utils import Units
from ..utils.ion_library import Ion

@dataclass
class TrappingParameters:
    """
    Dataclass to store parameters related to ion trapping.

    Parameters
    ----------
    ion : Optional[Ion], default=None
        The ion species and its properties.
    ion_height : Optional[float], default=None
        Height of the ion above the trap surface (in meters).
    rf_frequency : Optional[float], default=None
        RF frequency used for trapping (in Hz).
    rf_voltage : Optional[float], default=None
        RF voltage amplitude applied (in Volts).

    Examples
    --------
    >>> from ..utils.ion_library import Ion
    >>> params = TrappingParameters(ion=Ion('Yb', 171), ion_height=75e-6, rf_frequency=40e6, rf_voltage=150)
    """
    ion: Optional[Ion] = None
    ion_height: Optional[float] = None
    rf_frequency: Optional[float] = None
    rf_voltage: Optional[float] = None
    
class RadialDesigner:
    """
    Class for designing the radial geometry and operating parameters of a surface ion trap.

    Parameters
    ----------
    default_parameters : TrappingParameters
        Default set of trapping parameters for initialization.

    Examples
    --------
    >>> from py_iontrap_eda.designer.radial_designer import RadialDesigner, TrappingParameters
    >>> tp = TrappingParameters(ion=myIon, ion_height=75e-6, rf_frequency=40e6, rf_voltage=150)
    >>> rdes = RadialDesigner(tp)
    >>> design = rdes.compile(50e-6)
    >>> print(design.inner_dc_width, design.rf_width)
    """

    def __init__(self, default_parameters: TrappingParameters):
        """
        Initialize RadialDesigner.

        Parameters
        ----------
        default_parameters : TrappingParameters
            Default trapping parameters for this designer.
        """
        self.default_parameters = default_parameters
        self.default_rf_frequency = default_parameters.rf_frequency
        self.default_charge = default_parameters.ion.charge
        self.default_ion = default_parameters.ion
        self.default_mass = default_parameters.ion.ion_mass
        self.default_ion_height = default_parameters.ion_height

        self.depth_geometry_factor = lambda inner_dc_width, rf_width: (
            1 / np.pi**2 * rf_width**2
            / (2 * inner_dc_width + rf_width) ** 2
            / (
                np.sqrt(4 * inner_dc_width * (inner_dc_width + rf_width))
                + 2 * inner_dc_width
                + rf_width
            )
            ** 2
        )
        self.q_parameter_geometry_factor = lambda inner_dc_width, rf_width: (
            8 * rf_width
            / (2 * inner_dc_width + rf_width) ** 2
            / np.sqrt(inner_dc_width * (inner_dc_width + rf_width))
            / np.pi
        )
        self.secular_frequency_geometry_factor = lambda inner_dc_width, rf_width: (
            self.q_parameter_geometry_factor(inner_dc_width, rf_width) / (2 * np.sqrt(2))
        )
            
    def compile(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> RadialDesign:
        """
        Compile a RadialDesign instance from the specified geometry and trapping parameters.

        Parameters
        ----------
        inner_dc_width : float
            The width of the central DC electrode (in meters).
        trapping_parameters : Optional[TrappingParameters], default=None
            Trapping parameters to use. If None, uses the defaults from initialization.

        Returns
        -------
        RadialDesign
            An instance of RadialDesign with calculated parameters.

        Examples
        --------
        >>> design = rdes.compile(50e-6)
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width = self.get_rf_rail_width(inner_dc_width, trapping_parameters)
        return RadialDesign(inner_dc_width, rf_width)
        
    def compile_trapping_parameters(self, trapping_parameters: Optional[TrappingParameters] = None) -> TrappingParameters:
        """
        Fill in missing fields of trapping_parameters with defaults.

        Parameters
        ----------
        trapping_parameters : Optional[TrappingParameters], default=None
            The trapping parameters to fill.

        Returns
        -------
        TrappingParameters
            Updated trapping parameters with all fields populated.
        """
        if trapping_parameters is None:
            trapping_parameters = self.default_parameters
        if trapping_parameters.ion is None:
            trapping_parameters.ion = self.default_ion
        if trapping_parameters.ion_height is None:
            trapping_parameters.ion_height = self.default_ion_height
        if trapping_parameters.rf_frequency is None:
            trapping_parameters.rf_frequency = self.default_rf_frequency
        if trapping_parameters.rf_voltage is None:
            trapping_parameters.rf_voltage = self.default_rf_voltage
        return trapping_parameters
    
    def get_rf_rail_width(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> float:
        """
        Calculate the RF rail width given an inner DC width and trapping parameters.

        Parameters
        ----------
        inner_dc_width : float
            Inner DC electrode width (in meters).
        trapping_parameters : Optional[TrappingParameters], default=None
            Trapping parameters (optional).

        Returns
        -------
        float
            RF rail width (in meters).

        Examples
        --------
        >>> rf_width = rdes.get_rf_rail_width(50e-6)
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width = trapping_parameters.ion_height**2 / inner_dc_width - inner_dc_width
        return rf_width

    def trap_depth(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> float:
        """
        Calculate the trap depth for a given geometry and trapping parameters.

        Parameters
        ----------
        inner_dc_width : float
            Inner DC electrode width (in meters).
        trapping_parameters : Optional[TrappingParameters], default=None
            Trapping parameters (optional).

        Returns
        -------
        float
            Trap depth (in Joules).

        Examples
        --------
        >>> depth = rdes.trap_depth(50e-6)
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width = self.get_rf_rail_width(inner_dc_width, trapping_parameters)
        return (
            trapping_parameters.ion.charge**2 * trapping_parameters.rf_voltage**2 / (trapping_parameters.ion.ion_mass * trapping_parameters.rf_frequency**2) 
            * self.depth_geometry_factor(inner_dc_width, rf_width)
        )

    def q_parameter(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> float:
        """
        Calculate the (Mathieu) q-parameter for the given geometry and trapping parameters.

        Parameters
        ----------
        inner_dc_width : float
            Inner DC electrode width (in meters).
        trapping_parameters : Optional[TrappingParameters], default=None
            Trapping parameters (optional).

        Returns
        -------
        float
            q-parameter (dimensionless).

        Examples
        --------
        >>> q = rdes.q_parameter(50e-6)
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width = self.get_rf_rail_width(inner_dc_width, trapping_parameters)
        return (
            trapping_parameters.ion.charge * trapping_parameters.rf_voltage / (trapping_parameters.ion.ion_mass * trapping_parameters.rf_frequency**2) 
            * self.q_parameter_geometry_factor(inner_dc_width, rf_width)
        )

    def secular_frequency(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> float:
        """
        Calculate the secular oscillation frequency of the ion for the given geometry and trapping parameters.

        Parameters
        ----------
        inner_dc_width : float
            Inner DC electrode width (in meters).
        trapping_parameters : Optional[TrappingParameters], default=None
            Trapping parameters (optional).

        Returns
        -------
        float
            Secular frequency (in Hz).

        Examples
        --------
        >>> freq = rdes.secular_frequency(50e-6)
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width = self.get_rf_rail_width(inner_dc_width, trapping_parameters)
        return (
            trapping_parameters.ion.charge * trapping_parameters.rf_voltage / (trapping_parameters.ion.ion_mass * trapping_parameters.rf_frequency) * self.secular_frequency_geometry_factor(inner_dc_width, rf_width)
        )
    def geometry_sweep(
        self, 
        inner_dc_width_list: List[float], 
        trapping_parameters: Optional[TrappingParameters] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Sweep over a list of inner DC widths and calculate trap depth, q-parameter, and secular frequency.

        Parameters
        ----------
        inner_dc_width_list : list of float
            List of inner DC electrode widths (in meters) to sweep over.
        trapping_parameters : TrappingParameters, optional
            Trapping parameters to use for calculations. If None, uses default parameters.

        Returns
        -------
        rf_width_array : np.ndarray
            Array of calculated RF rail widths (in meters).
        trap_depth_array : np.ndarray
            Array of calculated trap depths (in eV).
        q_parameter_array : np.ndarray
            Array of calculated q-parameters (dimensionless).
        secular_frequency_array : np.ndarray
            Array of calculated secular frequencies (in Hz).

        Examples
        --------
        >>> inner_dc_widths = np.linspace(40e-6, 60e-6, 5)
        >>> rf_widths, trap_depths, q_params, sec_freqs = rdes.geometry_sweep(inner_dc_widths)
        >>> print(rf_widths)
        [0.123 0.125 0.127 0.129 0.130]
        >>> print(trap_depths)
        [0.123 0.125 0.127 0.129 0.130]
        >>> print(q_params)
        [0.123 0.125 0.127 0.129 0.130]
        >>> print(sec_freqs)
        [0.123 0.125 0.127 0.129 0.130]
        """
        trapping_parameters = self.compile_trapping_parameters(trapping_parameters)
        rf_width_array = np.array([
            self.get_rf_rail_width(inner_dc_width, trapping_parameters)
            for inner_dc_width in inner_dc_width_list
        ])
        trap_depth_array = np.array([
            self.trap_depth(inner_dc_width, trapping_parameters)
            for inner_dc_width in inner_dc_width_list
        ])
        q_parameter_array = np.array([
            self.q_parameter(inner_dc_width, trapping_parameters)
            for inner_dc_width in inner_dc_width_list
        ])
        secular_frequency_array = np.array([
            self.secular_frequency(inner_dc_width, trapping_parameters)
            for inner_dc_width in inner_dc_width_list
        ])
        return rf_width_array, trap_depth_array, q_parameter_array, secular_frequency_array
    def parameter_sweep(self, inner_dc_width: float, trapping_parameters: Optional[TrappingParameters] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Perform a parameter sweep for several trapping configurations.

        Parameters
        ----------
        inner_dc_width : float
            Inner DC electrode width (in meters).
        trapping_parameters : Optional[List[TrappingParameters]], default=None
            List of different trapping parameters to sweep over.

        Returns
        -------
        trap_depth_array : np.ndarray
            Array of trap depths.
        q_parameter_array : np.ndarray
            Array of q-parameters.
        secular_frequency_array : np.ndarray
            Array of secular frequencies.

        Examples
        --------
        >>> depth_array, q_array, freq_array = rdes.parameter_sweep(50e-6, [tp1, tp2, tp3])
        """
        trapping_parameters = [self.compile_trapping_parameters(tp) for tp in trapping_parameters]
        trap_depth_array = np.array([self.trap_depth(inner_dc_width, tp) for tp in trapping_parameters])
        q_parameter_array = np.array([self.q_parameter(inner_dc_width, tp) for tp in trapping_parameters])
        secular_frequency_array = np.array([self.secular_frequency(inner_dc_width, tp) for tp in trapping_parameters])
        return trap_depth_array, q_parameter_array, secular_frequency_array
