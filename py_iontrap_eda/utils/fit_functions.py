from typing import Callable, Any
from enum import Enum


class HarmonicFuncType(Enum):
    CENTER_POSITION = 1
    BASE_FORM = 2
    CENTERED = 3


class FitFunctions:
    @classmethod
    def harmonic_factory(cls, func_type: HarmonicFuncType):
        if func_type is HarmonicFuncType.CENTER_POSITION:

            def harmonic_center_position(x: float, a: float, x_0: float, y_0: float):
                return a * (x - x_0) ** 2 + y_0

            return harmonic_center_position

        elif func_type is HarmonicFuncType.BASE_FORM:

            def harmonic_base_form(x: float, a: float, b: float, c: float):
                return a * x**2 + b * x + c

            return harmonic_base_form

        elif func_type is HarmonicFuncType.CENTERED:

            def harmonic_centered(x: float, a: float, b: float):
                return a * x**2 + b

            return harmonic_centered
        else:
            raise ValueError(f"Unknown harmonic function type: {func_type}")
