"""Streamlit app smoke test: the script runs headless and raises nothing."""

from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")
pytest.importorskip("rocketcea")

APP = Path(__file__).resolve().parent.parent / "app.py"


def _app():
    return st_testing.AppTest.from_file(str(APP), default_timeout=300)


@pytest.mark.cea
def test_app_runs_with_defaults():
    at = _app().run()
    assert not at.exception
    assert not at.error
    assert len(at.tabs) == 7
    # The brief's example engine: delivered vacuum Isp 326.7 s.
    assert any("326.7 s" in m.value for m in at.metric)


@pytest.mark.cea
def test_app_reports_chart_limit_instead_of_extrapolating():
    at = _app().run()
    at.sidebar.radio[1].set_value("Fixed ε").run()
    eps_input = next(n for n in at.sidebar.number_input if n.label.startswith("ε"))
    eps_input.set_value(80.0).run()
    assert not at.exception
    assert at.error and "supply" in at.error[0].value.lower()


@pytest.mark.cea
def test_app_rp1_without_cp():
    at = _app().run()
    at.sidebar.selectbox[0].set_value("LOX/RP-1").run()
    assert not at.exception
    assert any(m.value == "n/a" for m in at.metric)
