"""Supplemental read-only numerical/temperature audits; no production mutations."""
import json
import os
from pathlib import Path
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(OUT / "mpl_cache"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(OUT / "numba_cache"))
import run_diagnostic as d
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def main():
    d.RESULTS.update(json.loads((OUT / "results.json").read_text(encoding="utf-8")))
    base = np.array(d.RESULTS["exact_corrected"]["G_s"])
    for temp in np.arange(123.5, 132.0001, .5):
        d.exact(f"exploratory_temperature_{temp:.2f}", dict(d.PARAMS, temp_c=float(temp)))
    off_temps = []
    for temp in np.arange(113., 123.0001, .5):
        tag = f"uncorrected_temperature_{temp:.2f}"
        d.exact(tag, dict(d.PARAMS, temp_c=float(temp)), gain_closure_enabled=False)
        off_temps.append(tag)
    eta_fast = d.displayed("ui_fast_eta0_5", detection_eff_pct=50.)
    noise_fast = d.displayed("ui_fast_noise_control", excess_noise=5., excess_noise_slope=1.)
    eta_balanced = d.displayed("ui_balanced_eta0_5", detection_eff_pct=50., resolution=d.fwm.FIDELITY_BALANCED)
    seed38 = d.exact("seed3_8_rounding_control", dict(d.PARAMS, probe_uw=3.8))
    compiled = d.exact("compiled_exact_corrected")
    # Expanded temperature range is an exploratory diagnostic, not a tolerance.
    def gain_at_temp(temp):
        cfg = d.kwargs(dict(d.PARAMS,temp_c=float(temp)))
        raw = d.fwm.compute_spectrum(d.PARAMS["opd"], scan_min=d.PARAMS["opd"]-3.039,
                                     scan_max=d.PARAMS["opd"]-3.030,coarse_points=10, **cfg)
        return raw["G_s"]
    match_last = brentq(lambda t: gain_at_temp(t)[-1]-d.MEAS[-1], 123., 132., xtol=1e-5)
    r_match = d.exact("temperature_matching_last_point", dict(d.PARAMS,temp_c=float(match_last)))
    optimum = minimize_scalar(lambda t: float(np.mean((gain_at_temp(t)-d.MEAS)**2)),
                              bounds=(113.,132.), method="bounded", options={"xatol":1e-4})
    d.exact("exploratory_best_temperature", dict(d.PARAMS,temp_c=float(optimum.x)))
    d.exact("peak_last_point_temperature",dict(d.PARAMS,temp_c=float(match_last)),
            eom_min=2.990,eom_max=3.050,points=241)
    audit = {}
    for tag,ref,readout in (("ui_fast_eta0_5","ui_fast",True),
                            ("ui_fast_noise_control","ui_fast",True),
                            ("ui_balanced_eta0_5","ui_balanced",True),
                            ("compiled_exact_corrected","exact_corrected",False),
                            ("seed3_8_rounding_control","exact_corrected",False)):
        a,b=d.RESULTS[tag],d.RESULTS[ref]
        ga,gb=np.array(a["G_s"]),np.array(b["G_s"])
        idx=int(np.argmax(np.abs(ga/gb-1)))
        entry=dict(max_relative_full_pct=float(100*np.max(np.abs(ga/gb-1))),
                   at_eom_GHz=a["eom_GHz"][idx],new_gain=ga[idx],reference_gain=gb[idx],
                   solved_nodes=a["metadata"]["response_estimator"]["solved_probe_points"],
                   reference_solved_nodes=b["metadata"]["response_estimator"]["solved_probe_points"])
        if readout:
            ga,gb=np.array(a["measured_eom_readout"]),np.array(b["measured_eom_readout"])
            entry["ten_point_max_relative_pct"]=float(100*np.max(np.abs(ga/gb-1)))
        audit[tag]=entry
    exploratory = {tag:d.fits(v["G_s"]) for tag,v in d.RESULTS.items()
                   if tag.startswith(("exploratory_temperature_","uncorrected_temperature_"))
                   or tag in ("temperature_matching_last_point","exploratory_best_temperature")}
    best_off=min(off_temps,key=lambda t:exploratory[t]["identity"]["rmse"])
    d.save_json("additional_checks.json",dict(audits=audit, exploratory_fits=exploratory,
                compiled_kernels_available=d.kernels.available(),
                temperature_matching_last_point_C=match_last,
                exploratory_best_common_temperature_C=float(optimum.x),
                best_uncorrected_temperature_in_5C=best_off,
                exploratory_temperature_is_not_tolerance=True))
    summary=json.loads((OUT/"summary.json").read_text(encoding="utf-8"))
    affine=summary["fits"]["exact_corrected"]["affine"]
    scale=summary["fits"]["exact_corrected"]["scale"]
    fig,ax=d.plt.subplots(1,2,figsize=(11,4.8),layout="constrained")
    ax[0].plot(d.EOM,d.MEAS,"ko-",label="Measured")
    ax[0].plot(d.EOM,scale["prediction"],"o-",label="Single output scale")
    ax[0].plot(d.EOM,affine["prediction"],"o-",label="Affine (diagnostic only)")
    ax[0].set_ylabel("Gain");ax[0].legend(fontsize=8)
    for y,label in ((d.MEAS,"Measured"),(base,"Corrected"),
                    (np.array(d.RESULTS["exact_uncorrected"]["G_s"]),"Uncorrected")):
        ax[1].plot(d.EOM,(y-y[0])/(y[-1]-y[0]),"o-",label=label)
    ax[1].set_ylabel("(Gain - first gain) / (last gain - first gain)")
    ax[1].legend(fontsize=8)
    for a in ax:
        a.set_xlabel("EOM frequency (GHz)");a.set_xlim(3.0395,3.0295);a.grid(alpha=.25)
    fig.savefig(OUT/"affine_diagnostic.png",dpi=180);d.plt.close(fig)
    print("Additional checks complete",json.dumps({"match_last_T":match_last,
          "best_T":float(optimum.x),"compiled":d.kernels.available()}),flush=True)


if __name__=="__main__":
    main()
