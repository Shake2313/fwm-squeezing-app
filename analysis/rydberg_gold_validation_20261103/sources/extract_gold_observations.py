"""Recover Gold figure coordinates with explicit calibration and source limits.

Fig.2 uses original PDF vector marker/path coordinates. Fig.5(b) is an embedded
raster; its marker centres and red dashed curve are recovered in pixel space.
No detector-voltage trace is relabelled as absolute optical transmission.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import fitz
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
PDF = REPO / "references/rydberg_experiments/sources/ju_2026_preprint.pdf"
POWERS = [.5, 1., 2., 4., 6., 8., 10.]


def write_json(name, value):
    (ROOT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2,
        allow_nan=False) + "\n", encoding="utf-8")


def write_csv(name, rows):
    with (ROOT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def linear_calibration(coordinates, values):
    coefficients = np.polyfit(coordinates, values, 1)
    residual = np.polyval(coefficients, coordinates) - values
    return {"coordinate": list(map(float, coordinates)), "value": list(map(float, values)),
            "map_value_equals_a_coordinate_plus_b": coefficients.tolist(),
            "max_tick_residual": float(np.max(np.abs(residual)))}


def calibrated(calibration, coordinate):
    return float(np.polyval(calibration["map_value_equals_a_coordinate_plus_b"], coordinate))


def path_points(drawing):
    points = []
    for item in drawing["items"]:
        if item[0] == "l":
            if not points:
                points.append(list(item[1]))
            points.append(list(item[2]))
    return np.array(points, dtype=float)


def source_width(x, y):
    order = np.argsort(x)
    x, y = x[order], y[order]
    baseline = float(np.median(y[np.abs(x) >= 2.2]))
    maximum = float(np.max(y[np.abs(x) < .3]))
    half = baseline + (maximum - baseline) / 2
    j = int(np.argmax(y))
    crossings = []
    for k in range(len(x) - 1):
        if (y[k] - half) * (y[k + 1] - half) <= 0 and y[k] != y[k + 1]:
            crossings.append(float(x[k] + (half-y[k])*(x[k+1]-x[k])/(y[k+1]-y[k])))
    left = max(value for value in crossings if value < x[j])
    right = min(value for value in crossings if value > x[j])
    return {"baseline_arb": baseline, "peak_arb": maximum,
            "peak_above_baseline_arb": maximum-baseline,
            "fwhm_mhz": right-left, "half_height_crossings_mhz": [left, right],
            "definition": "Median baseline at |frequency|>=2.2 MHz; peak maximum within 0.3 MHz; interpolated half-height crossings",
            "absolute_probe_power_uw": None,
            "absolute_optical_transmission": None}


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(PDF)
    page = doc[4]
    drawings = page.get_drawings()
    xcal = linear_calibration([378.74786377,396.66696167,414.58685303,
        432.50634766,450.42620850,468.34609985], [0,2,4,6,8,10])
    widthcal = linear_calibration([406.47100830,392.80444336,379.13784790,
        365.50115967,351.83459473,338.16799927], [1.5,2,2.5,3,3.5,4])
    peakcal = linear_calibration([322.25888062,304.04663086,285.83392334,
        267.62164307,249.41018677], [0,.2,.4,.6,.8])
    freqcal = linear_calibration([180.16879272,210.44288635,240.71603394,
        270.98919678,301.26239014], [-2,-1,0,1,2])
    signalcal = linear_calibration([381.00247192,356.73379517,332.46374512,
        308.19506836,283.92639160], [0,.1,.2,.3,.4])
    calibrations = {"fig2b_probe_power": xcal, "fig2b_width": widthcal,
                    "fig2b_peak": peakcal, "fig2a_frequency": freqcal,
                    "fig2a_signal": signalcal}
    points = []
    for quantity, ycal, yrange, unit in [("eit_fwhm", widthcal, (333,410), "MHz"),
                                       ("eit_peak", peakcal, (250,324), "arb.unit")]:
        for condition, color in [("compensated", (1.,0.,0.)), ("uncompensated", (0.,0.,1.))]:
            markers = []
            for index, drawing in enumerate(drawings):
                rect = drawing["rect"]
                if drawing["fill"] == color and 380 < rect.x0 < 472 and yrange[0] < rect.y0 < yrange[1] and rect.width < 3 and rect.height < 3:
                    markers.append((index, rect))
            markers.sort(key=lambda item: item[1].x0)
            if len(markers) != len(POWERS):
                raise RuntimeError((quantity, condition, len(markers)))
            for power, (index, rect) in zip(POWERS, markers):
                xc, yc = (rect.x0+rect.x1)/2, (rect.y0+rect.y1)/2
                points.append(dict(figure="2b", pdf_page=5, condition=condition,
                    probe_power_uw=power, digitized_probe_power_uw=calibrated(xcal, xc),
                    quantity=quantity, value=calibrated(ycal, yc), unit=unit,
                    pdf_x_pt=xc, pdf_y_pt=yc, drawing_index=index,
                    method="PDF vector marker rectangle centre",
                    comparison_resolution_bound=.03 if quantity == "eit_fwhm" else .005,
                    experimental_standard_uncertainty="not_assigned"))
    write_csv("fig2b_points.csv", points)
    curve_rows, shape_summary = [], {}
    for index in (48,49):
        drawing = drawings[index]
        condition = "compensated" if drawing["color"] == (1.,0.,0.) else "uncompensated"
        xy = path_points(drawing)
        x = np.array([calibrated(freqcal, value) for value in xy[:,0]])
        y = np.array([calibrated(signalcal, value) for value in xy[:,1]])
        shape_summary[condition] = source_width(x, y)
        for px, py, xx, yy in zip(xy[:,0],xy[:,1],x,y):
            curve_rows.append(dict(figure="2a", pdf_page=5, condition=condition,
                frequency_mhz=float(xx), transmission_signal_arb=float(yy),
                pdf_x_pt=float(px), pdf_y_pt=float(py), drawing_index=index))
    write_csv("fig2a_lineshape.csv", curve_rows)
    theory_rows = []
    for index in (299,300,371,372):
        drawing = drawings[index]
        condition = "compensated" if drawing["color"] == (1.,0.,0.) else "uncompensated"
        quantity = "eit_fwhm" if index < 301 else "eit_peak"
        ycal = widthcal if quantity == "eit_fwhm" else peakcal
        for px, py in path_points(drawing):
            theory_rows.append(dict(figure="2b", pdf_page=5, condition=condition,
                probe_power_uw=calibrated(xcal,px), quantity=quantity,
                value=calibrated(ycal,py), unit="MHz" if quantity == "eit_fwhm" else "arb.unit",
                drawing_index=index, method="PDF vector theoretical-fit polyline"))
    write_csv("fig2b_theory_curves.csv", theory_rows)

    # Recover raster with its PDF soft mask, composited onto white for QA.
    for pi, name in [(8,"figure4_source.png"),(9,"figure5_source.png")]:
        item = doc[pi].get_images()[0]
        pix = fitz.Pixmap(doc,item[0])
        if pix.alpha:
            pix = fitz.Pixmap(pix,0)
        if item[1]:
            pix = fitz.Pixmap(pix,fitz.Pixmap(doc,item[1]))
        image = Image.frombytes("RGBA",(pix.width,pix.height),pix.samples)
        background = Image.new("RGBA",image.size,"white")
        background.alpha_composite(image)
        background.convert("RGB").save(ROOT/name)
    image = np.asarray(Image.open(ROOT/"figure5_source.png").convert("RGB"))
    blue = (image[:,:,2]>150)&(image[:,:,0]<100)&(image[:,:,1]<130)
    red = (image[:,:,0]>180)&(image[:,:,1]<100)&(image[:,:,2]<100)
    # Major tick coordinates read from original embedded image, not page render.
    logxcal = linear_calibration([805.,1047.], [0.,1.])
    logycal = linear_calibration([375.,53.], [0.,3.])
    calibrations["fig5b_log10_probe_power_pixel"] = logxcal
    calibrations["fig5b_log10_sensitivity_pixel"] = logycal
    experimental = []
    for power in POWERS:
        xc = 805 + 242*np.log10(power)
        left, right = int(round(xc))-4, int(round(xc))+5
        profile = blue[130:330,left:right].sum(axis=1)
        rows = np.flatnonzero(profile == profile.max())+130
        yc = float(np.mean(rows))
        experimental.append(dict(figure="5b",pdf_page=10,probe_power_uw=power,
            sensitivity_nv_cm_sqrt_hz=10**calibrated(logycal,yc),
            image_x_px=float(xc),image_y_px=yc, method="Raster blue-marker densest-row centre; log-log axes",
            digitization_y_bound_px=2, experimental_standard_uncertainty="not_assigned"))
    write_csv("fig5b_experimental_points.csv",experimental)
    yy, xx = np.where(red & (np.indices(red.shape)[1]>=680) & (np.indices(red.shape)[0]>150))
    coefficients = np.polyfit(xx,yy,1)
    psn_curve = []
    for power in POWERS:
        xc = 805 + 242*np.log10(power)
        yc = float(np.polyval(coefficients,xc))
        psn_curve.append(dict(figure="5b",pdf_page=10,probe_power_uw=power,
            sensitivity_nv_cm_sqrt_hz=10**calibrated(logycal,yc), image_x_px=float(xc),
            image_y_px=yc, method="Fit red dashed raster centreline in log-log pixel space",
            digitization_y_bound_px=2, curve_kind="Published PSN theory, not new measurement"))
    write_csv("fig5b_psn_curve.csv",psn_curve)

    image4 = np.asarray(Image.open(ROOT/"figure4_source.png").convert("RGB"))
    gray = ((image4[:,:,0] > 40) & (image4[:,:,0] < 130)
            & (np.ptp(image4[:,:,:3],axis=2) < 4))
    figure4_points = []
    linear_y = linear_calibration([698.,424.], [0.,180.])
    calibrations["fig4d_field_pixel"] = linear_y
    for power in POWERS:
        xc = int(round(628 + 30*power))
        profile = gray[460:705,xc-3:xc+4].sum(axis=1)
        rows = np.flatnonzero(profile == profile.max())+460
        yc = float(np.mean(rows))
        figure4_points.append(dict(figure="4d", pdf_page=9, probe_power_uw=power,
            minimum_field_nv_cm=calibrated(linear_y,yc), image_x_px=xc,image_y_px=yc,
            method="Raster gray-marker densest-row centre; linear axes",
            digitization_y_bound_px=2, experimental_standard_uncertainty="not_assigned",
            unit_scope="Figure axis nV/cm; main text specifies 1 s integration and quotes field ASD"))
    write_csv("fig4d_experimental_points.csv",figure4_points)

    crop_specs = [(4,fitz.Rect(140,245,485,450),"fig2_source_crop.png")]
    for pi,name in [(8,"fig4_source_crop.png"),(9,"fig5_source_crop.png")]:
        rectangle = doc[pi].get_image_rects(doc[pi].get_images()[0][0])[0]
        rectangle = fitz.Rect(rectangle.x0-6,rectangle.y0-6,rectangle.x1+6,rectangle.y1+.5)
        crop_specs.append((pi,rectangle,name))
    for pi, rectangle, name in crop_specs:
        doc[pi].get_pixmap(matrix=fitz.Matrix(3,3),clip=rectangle).save(ROOT/name)
    protocols = [
        {"figure":"2a", "quantity":"Relative EIT spectral signal", "lo_rabi_mhz":0.,
         "lo_status":"EIT without AT dressing inferred from section ordering and single peak; explicit RF-off statement absent",
         "probe_power_uw":None,"coupling_power_mw":30.,"frequency_axis_identity":"NOT_REPORTED: probe scan versus coupling scan unspecified",
         "probe_detuning_definition":"Frequency relative to plotted EIT centre; absolute optical detuning absent",
         "absolute_transmission_available":False,"notes":"Unknown Pp: the 1.6 MHz spectrum cannot be assigned to 6 uW."},
        {"figure":"2b","quantity":"EIT peak and width vs Pp","lo_rabi_mhz":0.,
         "lo_status":"EIT measurement before LO/SIG section; no split peaks; explicit RF-off statement absent",
         "probe_power_uw":POWERS,"coupling_power_mw":30.,"coupling_rabi_mhz":3.,
         "frequency_axis_identity":"NOT_REPORTED: linewidth measured on an unidentified optical-frequency scan axis",
         "absolute_transmission_available":False,"notes":"Peaks are arbitrary signal units. Shape/calibration tests allowed; absolute transmission or normalized physical contrast not recoverable."},
        {"figure":"3","quantity":"AT spectrum and IF temporal signal","lo_rabi_mhz":3.7,
         "lo_carrier_ghz":37.07104,"if_khz":40.,"probe_power_uw":6.,
         "frequency_axis_identity":"NOT_REPORTED: probe scan versus coupling scan unspecified",
         "probe_detuning_definition":"Optimal spectral-slope operating point; absolute optical detuning not provided",
         "rf_detuning_status":"LO carrier reported but atomic zero reference not independently specified"},
        {"figure":"4","quantity":"Responsivity and measured sensitivity vs Pp","lo_rabi_mhz":3.7,
         "lo_status":"Main text supplies nominal LO 3.7 MHz and says LO optimized then fixed; another sentence says optimized conditions for each Pp. Missing Supplement Note 2 leaves per-Pp LO details unresolved.",
         "probe_power_uw":POWERS,"if_khz":40.,"probe_detuning_definition":"Operated near maximum spectral slope; numeric absolute optical point absent"},
        {"figure":"5b","quantity":"Measured and PSN sensitivity vs Pp","temperature_c":20.,
         "lo_rabi_mhz":3.7,"lo_status":"Inherited nominal heterodyne protocol; fixed/per-Pp LO details not fully specified",
         "probe_power_uw":POWERS,"if_khz":40.,"absolute_optical_detuning_mhz":None}]
    issues = [
        {"id":"linewidth_power_assignment","status":"SOURCE_CONDITION_UNRESOLVED",
         "finding":"Fig.2(a)/text width 1.6 MHz has no stated Pp; Fig.2(b) compensated width is about 1.60 MHz at 0.5 uW and 3.00 MHz at 6 uW. Do not use 1.6 MHz at 6 uW as a Gold validation anchor."},
        {"id":"absolute_transmission_and_contrast","status":"NOT_IDENTIFIABLE_FROM_MAIN_FIGURE",
         "finding":"Fig.2 axes are arbitrary signal units with a near-zero subtracted baseline; no dark/reference gain is supplied. Cannot infer absolute optical transmission, OD, or physical (Tpeak-Tbackground)/Tbackground contrast."},
        {"id":"psn_turnover_claim","status":"NARRATIVE_FIGURE_DISTINCTION",
         "finding":"The blue experimental curve has a minimum near 6 uW. The Fig.5(b) red PSN/SNL curve is monotonically decreasing approximately as Pp^-1/2; do not force a 6 uW PSN minimum from the narrative alone."},
        {"id":"LO_protocol","status":"SUPPLEMENT_REQUIRED_FOR_EXACT_PROTOCOL",
         "finding":"Nominal OmegaLO/2pi=3.7 MHz is reported. Main text discusses optimized conditions for each Pp and LO optimized then fixed. Exact per-Pp LO and numerical optical detunings remain unresolved without the Supplement."},
        {"id":"lineshape_width_definition", "status":"ESTIMATOR_AND_FIGURE_DIFFERENCE",
         "finding":"Direct half-height crossings of Fig.2(a) vector traces give red 1.625 MHz and blue 1.753 MHz, while the narrative quotes 1.6 and 1.9 MHz. No fit definition or Pp is supplied for this panel. Preserve both rather than treating the narrative as exact curve-derived FWHM."},
        {"id":"reported_vs_digitized_sensitivity", "status":"REPORTED_NUMBER_HAS_PRIORITY",
         "finding":"The narrative reports measured sensitivity 12.5(8) and PSN 11.2 nV/cm/sqrt(Hz) at 6 uW. Raster marker/line recovery gives about 13.4 and 11.8. Treat these as distinct approximate graphic extractions, with pixel uncertainty; do not overwrite reported numerical anchors."},
        {"id":"scan_axis_identity", "status":"NOT_REPORTED_FOR_GOLD",
         "finding":"The Gold main article labels the optical scan only Frequency (MHz), without identifying the scanned laser or absolute fixed detunings. Earlier Ju2025 thesis explicitly scans the probe with coupling fixed, but belongs to a different apparatus. Neither importing that protocol silently nor choosing a scan axis solely because its linewidth fits is justified."}]
    thesis_path = REPO / "analysis/rydberg_gold_map_20261103/sources/ju_2025_thesis.pdf"
    scan_audit = {"origin":"AI primary-source protocol audit",
        "gold_main":{"path":str(PDF.relative_to(REPO)).replace("\\","/"),
            "scan_laser":"not reported", "absolute_optical_detunings":"not reported",
            "locators":["PDF p.5 Fig.2 and caption", "PDF p.6 Fig.3 and caption"]},
        "prior_thesis":{"path":str(thesis_path.relative_to(REPO)).replace("\\","/"),
            "sha256":hashlib.sha256(thesis_path.read_bytes()).hexdigest(),
            "scan_laser":"probe", "coupling_frequency":"fixed on selected transition",
            "locators":["Printed p.16 / PDF p.25, 51D ladder", "Printed p.19 / PDF p.28, 40D ladder",
                        "Printed pp.17-18 / PDF pp.26-27, wavelength-axis scaling"],
            "applicability":"DIFFERENT_APPARATUS; supporting evidence only, does not establish Gold scan identity"},
        "recommendation":"Keep probe-scan and coupling-scan as separate conditional models until protocol is resolved; comparison success alone is not experimental identification."}
    write_json("scan_axis_audit.json",scan_audit)
    record = {"origin":"⟦AI⟧ Primary-figure extraction; not experimental recalibration",
        "source":{"path":str(PDF.relative_to(REPO)).replace("\\","/"),
                  "sha256":hashlib.sha256(PDF.read_bytes()).hexdigest(),"version":"arXiv:2606.04354v1","pages":15},
        "calibrations":calibrations,"protocols":protocols,"source_issues":issues,
        "fig2a_derived_shape":shape_summary,
        "fig2b_points":points,"fig4d_experimental_points":figure4_points,
        "fig5b_experimental_points":experimental,"fig5b_psn_curve":psn_curve,
        "reported_numeric_anchors":[
            {"quantity":"measured_field_asd","probe_power_uw":6.,"value":12.5,"printed_uncertainty":.8,
             "unit":"nV/cm/sqrt(Hz)","uncertainty_confidence":"not stated","locator":"PDF p.8 and abstract"},
            {"quantity":"PSN_field_asd","probe_power_uw":6.,"value":11.2,
             "unit":"nV/cm/sqrt(Hz)","locator":"PDF p.10"}],
        "fig5b_red_line_loglog_exponent":float(coefficients[0]*242*logycal["map_value_equals_a_coordinate_plus_b"][0]),
        "digitization_uncertainty_policy":"Vector comparison bounds 0.03 MHz/0.005 arb are conservative analyst comparison resolutions, not reported measurement standard uncertainties. Raster +/-2 pixel bounds describe digitization only. Published error-bar confidence is not assigned.",
        "source_figures_QA":"PDF pages 5,9,10 rendered and inspected; original marker/path and embedded-image extractions independently viewed",
        "supplement_status":"Main PDF references Notes 1 and 2; no ancillary link on official abstract. Both /src and /e-print return the identical main PDF. See supplement_search.json."}
    names=["fig2b_points.csv","fig2a_lineshape.csv","fig2b_theory_curves.csv",
           "fig4d_experimental_points.csv","fig5b_experimental_points.csv","fig5b_psn_curve.csv","extract_gold_observations.py",
           "fig2_source_crop.png","fig4_source_crop.png","fig5_source_crop.png",
           "scan_axis_audit.json","supplement_search.json"]
    record["artifact_sha256"]={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    write_json("gold_observations.json",record)
    print(json.dumps({"fig2b_points":len(points),"fig2a_shape":shape_summary,
        "fig5b_psn_exponent":record["fig5b_red_line_loglog_exponent"]},indent=2))


if __name__ == "__main__":
    main()
