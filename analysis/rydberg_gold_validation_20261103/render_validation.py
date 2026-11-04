"""Render saved Gold benchmark observations and hypotheses; no physical solves."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
FIELD = "sensitivity_nv_cm_sqrt_hz"


def read_json(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8"))


def read_csv(name):
    with (ROOT/name).open(encoding="utf-8-sig",newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key,value in row.items():
            try:
                row[key] = float(value)
            except (ValueError,TypeError):
                pass
    return rows


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    benchmark = read_json("benchmark_results.json")
    spectroscopy = read_csv("spectroscopy_comparison.csv")
    rf = read_csv("fixed_lo_rf.csv")
    observed = read_json("sources/gold_observations.json")
    gaussian = read_json("diagnostics/gaussian_trial.json")
    inputs = ["benchmark_results.json","spectroscopy_comparison.csv","fixed_lo_rf.csv",
              "sources/gold_observations.json","diagnostics/gaussian_trial.json",
              "sources/fig5b_experimental_points.csv","sources/fig5b_psn_curve.csv"]
    outputs = []
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
        "axes.titleweight":"bold","figure.facecolor":"white"})

    def save(figure,name):
        figure.savefig(OUT/name,dpi=180,bbox_inches="tight")
        plt.close(figure)
        outputs.append("figures/"+name)

    # A single model gain is shared between all powers and B-field conditions.
    fig,axes = plt.subplots(2,2,figsize=(11.8,8.8))
    fig.subplots_adjust(left=.075,right=.98,bottom=.20,top=.82,hspace=.28,wspace=.20)
    colors = {"off":"#7d7f87","on":"#246caf"}
    conditions = [("compensated","B-field compensated"),
                  ("uncompensated","B-field uncompensated")]
    train = set(benchmark["training_powers_uw"])
    for ci,(condition,label) in enumerate(conditions):
        source = sorted([row for row in spectroscopy if row["doppler"] == "on"
            and row["scan_axis"] == "Probe detuning" and row["condition"] == condition],
            key=lambda row:row["probe_power_uw"])
        if len(source) != 7:
            raise ValueError("Expected seven source powers per condition")
        for qi,(quantity,source_quantity,ylabel) in enumerate([
                ("eit_width_mhz","source_width_mhz","EIT FWHM (MHz)"),
                ("predicted_signal_arb","source_peak_arb","EIT peak signal (arb. unit)")]):
            axis = axes[ci,qi]
            for partition,marker,fill in [("training","o","#171717"),
                                           ("heldout","D","white")]:
                rows = [row for row in source if (row["probe_power_uw"] in train) == (partition == "training")]
                axis.scatter([row["probe_power_uw"] for row in rows],
                    [row[source_quantity] for row in rows],s=45,marker=marker,
                    facecolors=fill,edgecolors="#171717",linewidths=1.2,zorder=6,
                    label="Gold figure: "+partition)
            for mode in ["off","on"]:
                rows = sorted([row for row in spectroscopy if row["doppler"] == mode
                    and row["scan_axis"] == "Probe detuning" and row["condition"] == condition],
                    key=lambda row:row["probe_power_uw"])
                axis.plot([row["probe_power_uw"] for row in rows],[row[quantity] for row in rows],
                    color=colors[mode],lw=1.8,ls="--" if mode == "off" else "-",
                    label="CP uniform model: Doppler "+mode)
            if condition == "compensated":
                best = gaussian["selected_training_hypothesis"]
                rows = gaussian["selected_hypothesis_convergence_rows"]
                rows = sorted(rows,key=lambda row:row["probe_power_uw"])
                values = ([row["eit_fwhm_mhz"] for row in rows] if qi == 0 else
                    [row["relative_voltage_amplitude"]*best["one_global_amplitude_gain"] for row in rows])
                axis.plot([row["probe_power_uw"] for row in rows],values,color="#15916b",
                    ls="-.",lw=2.0,label="Gaussian HYPOTHESIS (rate bound)")
            else:
                axis.text(.03,.94,"Gaussian hypothesis not fitted to this condition",
                    transform=axis.transAxes,fontsize=8.5,va="top",color="#666666")
            axis.set(xlabel="Probe power (uW)",ylabel=ylabel,
                     title=label+" | "+("width" if qi == 0 else "relative signal"))
            axis.set_xticks([.5,2,4,6,8,10])
            axis.grid(alpha=.18)
            axis.set_ylim(bottom=0)
    handles,labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.085),
               ncol=3,frameon=False,fontsize=9)
    fig.suptitle("Gold EIT power-dependence benchmark\n"
        "LO off assumption | Pc = 30 mW | probe-scan models shown; Gold scan identity unresolved",
        fontsize=13,y=.96)
    fig.text(.5,.055,"Uniform models: one training-only signal gain shared across powers and B-field conditions. "
        "Gaussian: compensated training only.",ha="center",fontsize=8.7)
    fig.text(.5,.03,"Gaussian candidate selected at zero ground-rate bound; conditional hypothesis, not a validated model. "
        "No absolute-transmission calibration.",ha="center",fontsize=8.7,color="#9b2539")
    save(fig,"gold_eit_power_comparison.png")

    # Shape comparison uses its own six-microwatt normalization for every curve.
    fig,axis = plt.subplots(figsize=(9.7,6.3))
    fig.subplots_adjust(left=.095,right=.98,bottom=.24,top=.80)
    for policy,label,color in [("maximum_static_slope","CP thermal: maximum static slope","#784e9b"),
                               ("bare_two_photon_center","CP thermal: two-photon centre","#118f94")]:
        rows = sorted([row for row in rf if row["doppler"] == "on"
            and row["scan_axis"] == "Probe detuning" and row["detuning_policy"] == policy],
            key=lambda row:row["probe_power_uw"])
        denominator = next(row[FIELD] for row in rows if row["probe_power_uw"] == 6.)
        axis.semilogy([row["probe_power_uw"] for row in rows],
            [row[FIELD]/denominator for row in rows],".-",color=color,label=label,lw=1.8,ms=7)
    for name,label,color,style in [
            ("sources/fig5b_experimental_points.csv","Gold measured shape (figure-derived)","#202020","o-"),
            ("sources/fig5b_psn_curve.csv","Gold plotted PSN shape (figure-derived)","#d43d45","--")]:
        rows = sorted(read_csv(name),key=lambda row:row["probe_power_uw"])
        denominator = next(row[FIELD] for row in rows if row["probe_power_uw"] == 6.)
        axis.semilogy([row["probe_power_uw"] for row in rows],
            [row[FIELD]/denominator for row in rows],style,color=color,label=label,lw=1.7,ms=5)
    axis.axhline(1,color=".65",lw=.8)
    axis.axvline(6,color=".65",lw=.8,ls=":")
    axis.set(xlabel="Probe power (uW)",ylabel="Field sensitivity / own value at 6 uW",
             title="Fixed LO Rabi / 2pi = 3.7 MHz | IF = 40 kHz | T = 20 C | Pc = 30 mW")
    axis.set_xticks([.5,1,2,4,6,8,10])
    axis.grid(alpha=.2,which="both")
    axis.legend(fontsize=9,loc="upper center",framealpha=.95)
    fig.suptitle("Probe-power response at fixed LO: normalized shapes\n"
        "Different declared optical-point policies; exact Gold optical detunings are unreported",
        fontsize=13,y=.96)
    fig.text(.5,.12,"All curves normalized independently at 6 uW. "
        "Graphic extraction is not substituted for the paper's reported absolute anchors:",ha="center",fontsize=9)
    fig.text(.5,.085,"measured 12.5(8), PSN 11.2 nV/cm/sqrt(Hz). "
        "The plotted red PSN curve decreases with power.",ha="center",fontsize=9)
    fig.text(.5,.045,"CP curves assume zero technical noise. This tests shape under declared policies; "
        "it does not establish matched absolute operating detunings.",ha="center",fontsize=9,color="#9b2539")
    save(fig,"fixed_lo_normalized_power_response.png")

    guard_files = [name for name in ("current_resolution_audit.json","refined_resolution_audit.json",
                                    "refined_current_resolution_audit.json","refined_current_resolution.json")
                   if (ROOT/name).exists()]
    if guard_files:
        # Each file must declare its own new guards; unlike repeated training
        # points these maxima test the numerical interpolation policy.
        records = []
        for name in guard_files:
            document = read_json(name)
            for row in document.get("reports",[]):
                row = dict(row)
                if name == "refined_resolution_audit.json":
                    row["density_normalized"] = True
                records.append((name,row))
        if records:
            inputs.extend(guard_files)
            fig,axes = plt.subplots(1,2,figsize=(11.5,5.6))
            fig.subplots_adjust(left=.07,right=.98,bottom=.27,top=.76,wspace=.22)
            for axis,policy,title in zip(axes,["maximum_static_slope","bare_two_photon_center"],
                                        ["Maximum static slope","Two-photon centre"]):
                rows = [(name,row) for name,row in records if row["detuning_policy"] == policy]
                heights = [100*row["maximum_relative_error"] for _,row in rows]
                labels = [("V1\nraw log S" if not row.get("density_normalized",False) else
                           "V1\ndensity residual") if name == "current_resolution_audit.json" else
                          "V2\ndensity + Pc nodes" for name,row in rows]
                bars = axis.bar(range(len(rows)),heights,color=["#5c8bba" if row.get("density_normalized",False)
                    else "#9b9ea6" for _,row in rows],width=.62)
                for bar,height in zip(bars,heights):
                    inside = .7 < height < 1.
                    axis.text(bar.get_x()+bar.get_width()/2,height-.04 if inside else height+.035,
                        f"{height:.3f}%",ha="center",va="top" if inside else "bottom",
                        fontsize=9,color="white" if inside else "black")
                axis.axhline(1,color="#b32c43",ls="--",lw=1.2,label="1% numerical criterion")
                axis.set_xticks(range(len(rows)),labels,fontsize=9)
                axis.set(ylabel="Maximum absolute guard error (%)",title=title,
                         ylim=(0,max([1.3]+heights)*1.22))
                axis.grid(axis="y",alpha=.15)
                axis.legend(fontsize=8.5,frameon=False)
            fig.suptitle("Sparse-prediction guard checks in the corrected CP model\n"
                "Known density factor separated from the residual response",fontsize=13,y=.96)
            fig.text(.5,.13,"V1 and V2 use separately predeclared guard sets. "
                "A failed V1 policy remains a recorded result.",ha="center",fontsize=9)
            fig.text(.5,.08,"Passing a finite guard set supports a bounded numerical predictor. "
                "It is not global convergence or experimental validation.",ha="center",fontsize=9,color="#9b2539")
            save(fig,"sparse_prediction_guard_errors.png")
    manifest = {"origin":"AI rendering of saved Gold validation data",
        "scope":"No new physical solves; no implied experimental calibration PASS",
        "input_sha256":{name:digest(ROOT/name) for name in inputs},
        "code_sha256":digest(Path(__file__)),
        "output_sha256":{name:digest(ROOT/name) for name in outputs}}
    (ROOT/"figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"figures":outputs},indent=2))


if __name__ == "__main__":
    main()
