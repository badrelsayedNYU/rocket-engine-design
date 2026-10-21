"""Stage 6: engine data integrity and validation pipeline."""
from pathlib import Path
import sys

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import run_validation as rv  # noqa: E402

from engine_design.losses import bell_divergence_factor  # noqa: E402


def test_every_cited_source_exists_and_values_are_labeled():
    raw = yaml.safe_load(rv.DATA.read_text())
    for engine in raw["engines"]:
        for key in ("p_c", "mixture_ratio", "area_ratio", "isp_vac", "thrust_vac"):
            entry = engine[key]
            if entry.get("value") is None:
                continue
            assert entry["quality"] in ("primary", "secondary", "derived", "assumption")
            for source in entry.get("source", []):
                assert source in raw["sources"], f"{engine['name']} {key}: {source}"


def test_unit_conversion_to_si():
    engines = {e.name: e for e in rv.load_engines()}
    # NIST SP 811: 2994 psi = 20.643 MPa; 24,750 lbf = 110.093 kN.
    assert engines["RS-25"].p_c_Pa == pytest.approx(2994 * 6894.757293168)
    assert engines["RS-25"].p_c_Pa / 1e6 == pytest.approx(20.643, abs=1e-3)
    assert engines["RL10B-2"].thrust_vac_N / 1e3 == pytest.approx(110.093, abs=1e-3)
    assert engines["Merlin 1D"].p_c_Pa == 9.7e6


def test_derived_sea_level_isp():
    # Equal mass flow at both conditions: Isp_sl = Isp_vac * F_sl / F_vac.
    rs25 = {e.name: e for e in rv.load_engines()}["RS-25"]
    assert rs25.isp_sl_s == pytest.approx(452.3 * 418000 / 512300)


def test_charts_never_extrapolated():
    engines = {e.name: e for e in rv.load_engines()}
    for name in ("RS-25", "RL10B-2"):
        with pytest.raises(ValueError):
            bell_divergence_factor(engines[name].eps, 0.8)
        assert rv.divergence(engines[name])[1].startswith("assumption")
    assert rv.divergence(engines["Merlin 1D"])[0] == bell_divergence_factor(16, 0.8)


@pytest.fixture(scope="module")
def predictions():
    pytest.importorskip("rocketcea")
    return {e.name: rv.predict(e) for e in rv.load_engines()}


@pytest.mark.cea
def test_ideal_shifting_exceeds_published(predictions):
    # docs/spec/03_VALIDATION.md: below published means something is wrong.
    for p in predictions.values():
        assert p.shifting.Isp_vac_s > p.engine.isp_vac_s
        assert p.frozen.Isp_vac_s < p.shifting.Isp_vac_s
        assert p.delivered.Isp_vac_s < p.shifting.Isp_vac_s


@pytest.mark.cea
def test_validation_path_matches_design_engine(predictions):
    # The validation route is the same physics as design_engine.
    from engine_design import EngineRequirements, design_engine
    p = predictions["Merlin 1D"]
    e = p.engine
    d = design_engine(EngineRequirements(propellants=e.propellants, thrust_N=845e3,
                                         thrust_condition="sea_level", p_c_Pa=e.p_c_Pa,
                                         mixture_ratio=e.OF, nozzle=e.eps))
    assert d.delivered.Isp_vac_s == pytest.approx(p.delivered.Isp_vac_s, rel=1e-12)
    assert d.ideal_frozen.Isp_vac_s == pytest.approx(p.frozen.Isp_vac_s, rel=1e-12)


@pytest.mark.cea
def test_eta_c_star_scales_delivered_isp_exactly(predictions):
    e = predictions["RL10B-2"].engine
    base = rv.predict(e).delivered.Isp_vac_s
    assert rv.predict(e, eta_c_star=0.98 * 1.01).delivered.Isp_vac_s == pytest.approx(base * 1.01)


@pytest.mark.cea
def test_merlin_is_illustrative_band(predictions):
    p = predictions["Merlin 1D"]
    assert p.engine.illustrative and p.engine.OF_band == (2.3, 2.4)
    assert "Merlin 1D" not in rv.TARGET_PERCENT
    lo, hi = p.isp_range("delivered")
    assert lo < hi and lo <= p.delivered.Isp_vac_s <= hi
    for of in (2.3, 2.4):
        isp = rv.predict(rv.replace(p.engine, OF=of, OF_band=None)).delivered.Isp_vac_s
        assert lo <= isp <= hi
    assert "illustrative" in rv._verdict(p)
    # Non-illustrative engines have no band.
    rs25 = predictions["RS-25"]
    assert rs25.isp_range("delivered") == (rs25.delivered.Isp_vac_s,) * 2


@pytest.mark.cea
def test_markdown_table(predictions):
    text = rv.markdown(list(predictions.values()), {})
    assert "| Engine | Cycle | Published | Predicted (frozen) | Predicted (shifting) " in text
    for name in predictions:
        assert f"| {name} |" in text
