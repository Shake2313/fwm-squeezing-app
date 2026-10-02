"""Cross-document claims that must change together."""
import json
from pathlib import Path
from urllib.parse import quote


ROOT = Path(__file__).resolve().parent.parent


def _read(relative_path):
    return (ROOT / relative_path).read_text(encoding="utf-8")


def _flat(relative_path):
    return " ".join(_read(relative_path).split())


def test_repository_visibility_wording_is_consistently_public():
    readme = _read("README.md")
    guide = _read("docs/Userguide/GABES_User_Guide_v2.html")
    assert "github.com/Shake2313/fwm-squeezing-app` (public)" in readme
    assert "로컬 설치용 공개 소스" in guide
    assert "fwm-squeezing-app` (private)" not in readme
    assert "소스 (private)" not in guide


def test_archived_paper_preserves_reproducibility_and_zeeman_caveats():
    paper_flat = _flat("docs/squeezing_report/squeezing_paper.tex")
    assert "openly available in the GitHub repository" in paper_flat
    assert "retained raw scan archives" in paper_flat
    assert "available only in their archived rendered form" in paper_flat
    assert "archived scan data underlying every figure" not in paper_flat
    assert "scripts that regenerate every figure" not in paper_flat
    assert paper_flat.count("24-level Zeeman CG-sum diagnostic") >= 2
    assert "CG-weighted Zeeman corrections" not in paper_flat


def test_readme_preserves_biphoton_compatibility_and_coupling_limits():
    # Generic biphoton/API claims belong to the active project documentation;
    # the current analytic publication covers seeded double-Lambda theory.
    readme_flat = _flat("README.md")
    assert "Cs cascade channels land within ~30 %" not in readme_flat
    assert "OD/cell length do not directly scale the rate" in readme_flat
    assert "intrinsic source FWHM" in readme_flat
    assert (
        "legacy `fwhm_ns` API key is retained as an alias for the detected width"
    ) in readme_flat
    assert "source FWHM is about 0.17 ns" in readme_flat
    assert "detected FWHM to about 0.50 ns" in readme_flat
    assert "pair-rate scale remains reference-anchored" in readme_flat
    assert "That agreement is not an independent validation" in readme_flat
    assert "reference residual × additional mode-overlap penalty" in readme_flat
    assert "no unsupported numerical split of `0.74`" in readme_flat
    assert "separate from Ultra's normalized axial crossing-angle profile" in (
        readme_flat
    )


def test_current_publications_preserve_conditional_claim_boundaries():
    # Publication policy makes the manifest authoritative; retired editions
    # must not become prerequisites for checking the current scientific claims.
    publications = json.loads(
        _read("docs/grand_challenge/current_publications.json"))["publications"]
    readme = _read("README.md")
    for publication in publications.values():
        assert (ROOT / publication["source"]).is_file()
        assert (ROOT / publication["artifact"]).is_file()
        assert f'({quote(publication["artifact"])})' in readme

    physics = _flat(publications["analytic_reconstruction"]["source"])
    report = _flat(publications["squeezing_report"]["source"])
    quotient = _flat(publications["quotient_structure"]["artifact"])
    assert "$0.74$ remains an assumed residual" in physics
    assert "not a microscopic squeezing prediction" in physics
    assert "not certificates of an experimental squeezing value" in physics
    assert "The reported thermal ensemble is not converged" in physics
    assert "원자 입력의 독립 예측 아님" in report
    assert "전체 thermal 수렴 미인증" in report
    assert "가정, 정확한 동치, 정보의 손실, 수치 근사와 실험 검증을 구별" in quotient
    assert "실제 장치 예측 미인증" in quotient
