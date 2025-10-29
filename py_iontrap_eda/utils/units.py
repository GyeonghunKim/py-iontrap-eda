from typing import Optional


class Units:
    # SI units
    m = 1
    kg = 1
    s = 1
    A = 1

    # Length units
    cm = 1e-2 * m
    mm = 1e-3 * m
    um = 1e-6 * m
    nm = 1e-9 * m
    inch = 2.54 * cm
    mil = 1e-3 * inch

    # Frequency units
    Hz = 1 / s
    kHz = 1e3 * Hz
    MHz = 1e6 * Hz
    GHz = 1e9 * Hz

    # Charge units
    C = 1

    # Energy units
    J = 1

    # Potential units
    V = J / C
    mV = 1e-3 * V
    eV: Optional[float] = None
    meV: Optional[float] = None


class Constants:
    amu = 1.66053906660e-27 * Units.kg
    e = 1.602176634e-19 * Units.C


Units.eV = Constants.e * Units.V
Units.meV = 1e-3 * Units.eV
