"""
Object classes for labels used to for supervised-training
"""

from typing import Tuple


class Label:
    def __init__(self, value, functional: str, basis: str):
        self.value = value
        self.functional = functional
        self.basis = basis
        return


class Potential(Label):
    def __init__(self, value: float, functional: str, basis: str):
        super().__init__(value, functional, basis)

    def __repr__(self):
        return f"({self.value:.7f}, {self.functional}, {self.basis})"

    def export_json(self):
        key = ".".join(("potential", self.functional, self.basis))
        return key, self.value


class Forces(Label):
    def __init__(self, value: Tuple[float, float, float], functional: str, basis: str):
        super().__init__(value, functional, basis)

    def __repr__(self):
        return f"([{self.value[0]:.7f}, {self.value[1]:.7f}, {self.value[2]:.7f}], {self.functional}, {self.basis})"

    def export_json(self):
        key = ".".join(("forces", self.functional, self.basis))
        return key, self.value


class Dipole(Label):
    def __init__(self, value: Tuple[float, float, float], functional: str, basis: str):
        super().__init__(value, functional, basis)

    def __repr__(self):
        return f"([{self.value[0]:.7f}, {self.value[1]:.7f}, {self.value[2]:.7f}], {self.functional}, {self.basis})"

    def export_json(self):
        key = ".".join(("dipole", self.functional, self.basis))
        return key, self.value


class Quadrupole(Label):
    def __init__(self, value: Tuple[float, float, float, float, float, float], functional: str, basis: str):
        super().__init__(value, functional, basis)

    def __repr__(self):
        return f"([{self.value[0]:.7f}, {self.value[1]:.7f}, {self.value[2]:.7f}, {self.value[3]:.7f}," \
               f" {self.value[4]:.7f}, {self.value[5]:.7f}], {self.functional}, {self.basis})"

    def export_json(self):
        key = ".".join(("quadrupole", self.functional, self.basis))
        return key, self.value


class Charge(Label):
    def __init__(self, value: float, partitioning: str, functional: str, basis: str):
        super().__init__(value, functional, basis)
        self.partitioning = partitioning

    def __repr__(self):
        return f"({self.value:.4f}, {self.partitioning}, {self.functional}, {self.basis})"

    def export_json(self):
        key = ".".join(("charge", self.partitioning, self.functional, self.basis))
        return key, self.value
