"""Cross-document claims that must change together."""
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def test_repository_visibility_wording_is_consistently_public():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    guide = (
        ROOT / "docs" / "Userguide" / "GABES_User_Guide_v2.html"
    ).read_text(encoding="utf-8")
    paper = (ROOT / "docs" / "squeezing_report" / "squeezing_paper.tex").read_text(
        encoding="utf-8")
    physics = (
        ROOT
        / "docs"
        / "FWM physics and analytic reconstruction"
        / "FWM_physics.tex"
    ).read_text(encoding="utf-8")
    readme_flat = " ".join(readme.split())
    paper_flat = " ".join(paper.split())
    physics_flat = " ".join(physics.split())

    assert "github.com/Shake2313/fwm-squeezing-app` (public)" in readme
    assert "로컬 설치용 공개 소스" in guide
    assert "openly available in the GitHub repository" in paper_flat
    assert "retained raw scan archives" in paper_flat
    assert "available only in their archived rendered form" in paper_flat
    assert "archived scan data underlying every figure" not in paper_flat
    assert "scripts that regenerate every figure" not in paper_flat
    assert paper_flat.count("24-level Zeeman CG-sum diagnostic") >= 2
    assert "CG-weighted Zeeman corrections" not in paper_flat
    assert "Cs cascade channels land within ~30 %" not in readme_flat
    assert "OD/cell length do not directly scale the rate" in readme_flat
    assert "intrinsic source FWHM" in readme_flat
    assert "legacy `fwhm_ns` API key" in readme_flat
    assert r"source} width converges to $\approx0.17$" in physics_flat
    assert r"detected} width of $\approx0.50$" in physics_flat
    assert r"legacy \code{fwhm\_ns} key remains an alias" in physics_flat
    assert "not an exact dual-observable calibration" in physics_flat
    assert "reference residual × additional mode-overlap penalty" in readme_flat
    assert "no unsupported numerical split of `0.74`" in readme_flat
    assert "does not invent one" in physics_flat
    assert "separate from Ultra's normalized axial crossing-angle profile" in (
        physics_flat
    )
    assert "fwm-squeezing-app` (private)" not in readme
    assert "소스 (private)" not in guide
