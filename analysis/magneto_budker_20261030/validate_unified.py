"""Verify the actual registered scheme, scan parity and former curve gaps."""
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gabes import schemes
from gabes.schemes.magneto_unified import FIELD_SCAN, OPTICAL_SCAN
from gabes_ui.export import to_json_primitive


def run():
    names = ["gabes/magneto_multilevel.py", "gabes/schemes/magneto.py",
             "gabes/schemes/magneto_unified.py", "gabes/magneto_noise.py"]
    snapshot = {n: hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}
    target = Path(__file__).parent / "results"
    target.mkdir(parents=True, exist_ok=True)
    scheme = schemes.get("magneto")
    base = dict(scheme.defaults(), show_shot_noise=False)
    settings = {
        "default_field": base,
        "detuned_field": dict(base, laser_detuning_mhz=-300.),
        "optical": dict(base, scan_axis=OPTICAL_SCAN),
        "optical_refined": dict(base, scan_axis=OPTICAL_SCAN, velocity_points=10,
                                detuning_points=97),
    }
    evidence, raw_cases = {}, {}
    for label, params in settings.items():
        started = time.perf_counter()
        raw = scheme.compute(params)
        view = scheme.observables(raw, params)
        evidence[label] = dict(params=params, seconds=time.perf_counter()-started,
                                readout=view["magneto_readout"])
        raw_cases[label] = raw
        view["figure"].savefig(target/f"unified_{label}.png", dpi=155)
        plt.close(view["figure"])
        print(label, f"{evidence[label]['seconds']:.2f}s", flush=True)
    b = evidence["detuned_field"]["readout"]
    optical = evidence["optical"]["readout"]
    refined = evidence["optical_refined"]["readout"]
    bi = b["operating_index"]
    oi = optical["x"].index(-300.)
    np.testing.assert_allclose(b["rotation_slope_rad_per_t"][bi],
                                optical["rotation_slope_rad_per_t"][oi], rtol=1e-10)
    np.testing.assert_allclose(b["transmission"][bi], optical["transmission"][oi], rtol=1e-12)
    gaps = [optical["x"].index(value) for value in (0., 300., 400.)]
    assert all(optical["slope_valid"][i] for i in gaps)
    assert all(optical["sensitivity"]["field_asd_t_per_sqrt_hz"][i] is not None for i in gaps)
    assert all(optical["slope_valid"])
    assert all(v is not None for v in optical["sensitivity"]["field_asd_t_per_sqrt_hz"])
    schema = {s.name for s in scheme.param_schema()}
    assert "calculation" not in schema and not any(n.startswith("budker_") for n in schema)
    assert base["scan_axis"] == FIELD_SCAN
    assert all(type(schemes.get(name)) is type(scheme) for name in ("hanle", "eia", "nmor"))
    delta = np.asarray(optical["x"])
    indices = np.array([refined["x"].index(value) for value in delta])
    slopes = np.asarray(optical["rotation_slope_rad_per_t"])
    fine_slopes = np.asarray(refined["rotation_slope_rad_per_t"])[indices]
    assessment = dict(
        all_registered_aliases_same_multihyperfine_engine=True,
        default_B_scan_uses_Budker_engine=True,
        no_model_selection_control=True,
        scan_axes_agree_at_same_physical_conditions=True,
        former_gap_detunings_mhz=[0., 300., 400.],
        former_gaps_all_have_finite_valid_analytic_slope_and_asd=True,
        optical_finite_asd_samples=len(delta),
        slope_max_relative_residual=max(optical["slope_relative_residual"]),
        shared_delta_slope_error_over_peak=float(np.max(np.abs(slopes-fine_slopes))/np.max(np.abs(fine_slopes))),
        optimum_detuning_mhz=refined["optimum_detuning_mhz"],
        optimum_asd_t_per_sqrt_hz=refined["optimum_field_asd_t_per_sqrt_hz"],
        optimum_relative_change=(refined["optimum_field_asd_t_per_sqrt_hz"]/
                                 optical["optimum_field_asd_t_per_sqrt_hz"]-1),
        status="conditional effective model; not quantitative paper reproduction",
    )
    assert snapshot == {n: hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in names}, "source changed during validation"
    output = dict(source_sha256=snapshot,
                  cases=evidence, assessment=assessment)
    (target/"unified_validation.json").write_text(
        json.dumps(to_json_primitive(output), ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8")
    print(json.dumps(assessment, ensure_ascii=False, indent=2), flush=True)
    return output


if __name__ == "__main__":
    run()
