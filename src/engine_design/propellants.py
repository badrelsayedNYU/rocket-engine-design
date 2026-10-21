"""Propellant labels mapped to built-in RocketCEA reactant cards.

Source: https://rocketcea.readthedocs.io/en/latest/propellants.html
The cards include reactant state and enthalpy, not just chemical composition.
Acceptance by the installed backend is checked in tests/test_propellants.py.
"""

from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class PropellantPair:
    """Display name and RocketCEA card identifiers (no performance assumptions)."""

    name: str
    cea_oxidizer: str
    cea_fuel: str


PROPELLANTS = MappingProxyType({
    "LOX/LH2": PropellantPair("LOX/LH2", "LOX", "LH2"),
    "LOX/RP-1": PropellantPair("LOX/RP-1", "LOX", "RP1"),
    "LOX/CH4": PropellantPair("LOX/CH4", "LOX", "CH4"),
})
