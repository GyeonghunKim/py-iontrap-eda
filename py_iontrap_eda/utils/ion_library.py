from .units import Units, Constants


class Ion:
    def __init__(
        self, name: str, atomic_number: int, mass_number: int, relative_mass: float
    ):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.relative_mass = relative_mass
        self.ion_mass = relative_mass * Constants.amu
        self.charge = Constants.e


class IonLibrary:
    Ba138 = Ion("Ba138", 56, 138, 137.905247)
    Yb171 = Ion("Yb171", 70, 171, 170.9363258)
    Ba137 = Ion("Ba137", 56, 137, 136.905827)
    Ca40 = Ion("Ca40", 20, 40, 39.9625909)