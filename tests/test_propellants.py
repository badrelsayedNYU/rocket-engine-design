"""Stage 0 integration check: names must work in a real CEA calculation."""

import math

import pytest

from engine_design.propellants import PROPELLANTS


@pytest.mark.cea
@pytest.mark.parametrize("name,of_ratio", [
    ("LOX/LH2", 6.0),
    ("LOX/RP-1", 2.5),
    ("LOX/CH4", 3.5),
])
def test_cea_accepts_propellant_pair(name, of_ratio, tmp_path, monkeypatch):
    """Exercise card parsing and Fortran; these inputs are not engine data.

    A positive finite output is an integration check, not physics validation.
    The units wrapper makes chamber pressure explicit in SI.
    Reference: https://rocketcea.readthedocs.io/en/latest/simple_examples.html
    """
    pytest.importorskip("rocketcea")
    from rocketcea import cea_obj
    from rocketcea.cea_obj_w_units import CEA_Obj

    # RocketCEA 1.2.3 exposes this directory as a module global.
    # Keep Fortran scratch files out of the user's home directory.
    monkeypatch.setattr(cea_obj, "ROCKETCEA_DATA_DIR", str(tmp_path))
    pair = PROPELLANTS[name]
    cea = CEA_Obj(
        oxName=pair.cea_oxidizer,
        fuelName=pair.cea_fuel,
        pressure_units="Pa",
        isp_units="sec",
    )
    isp_s = cea.get_Isp(Pc=10e6, MR=of_ratio, eps=40.0)
    assert math.isfinite(isp_s) and isp_s > 0
