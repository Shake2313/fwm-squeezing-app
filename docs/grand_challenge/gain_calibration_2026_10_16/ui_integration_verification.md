# ⟦AI⟧ Measured-reference UI integration verification

Verification notes, 2026-10-16. This records implementation checks; it does not add a user goal or claim independent physical validation.

- The conditional lab preset adopts the external incident seed proxy explicitly, preserves its full numeric precision, and leaves the measured-reference selector off. Enabling it requires the recorded optical context and RF support. Detector efficiency and classical seed noise remain readout controls.
- Current measured gains use the shared deterministic log-PCHIP helper. They occupy a separate `measured_gain_reference` ledger and reference-labelled point metrics. Atomic gains, canonical transfer and NF remain separate. The NF caption explicitly states that the measured gain interpolation does not recalibrate squeezing.
- The appended comparison view contains 23 raw points and their interpolation for each arm. The existing full-scan view remains index 0 with its original identifier. [Rendered comparison](measured_gain_ui_comparison.png) was inspected: labels and both arms are legible, with node agreement identified as construction rather than independent validation.
- Machine roundoff reconstructed at RF endpoints is normalized only within 16 floating-point spacings of 3022 or 3044 MHz. All 23 raw-axis roundtrips retain support and both-arm node agreement; a 0.00000001 MHz outward offset remains unavailable. A plus-branch operating-point request cannot reuse the lower-seed reference.
- Strict JSON export retains profile, RF domain, status, support and both gain arrays, including unavailable values represented without nonstandard NaN literals.

Targeted verification: 52 passed for the new integration file plus `tests/test_fwm_generic.py`; 157 passed for the integration, excess-noise, passive-L1 integration, historical-gain, TPD-gain, fast-tier and scheme-render files. The two warnings in the broader run concern the existing UTC timestamp export API.

The Streamlit AppTest exercised the real preset callback, new checkbox, rail controls, derived values, metrics and operating-point table across reruns. It confirmed the initial off state, unavailable unmatched context, preset reset to off, eligible opt-in, and the exact 3.921545933233986 µW proxy surviving the power slider. The preexisting formatted segmented mode control was omitted only from this harness because Streamlit AppTest cannot serialize it on rerun; the production mode control is unchanged.

Python changes for this integration were confined to `gabes/schemes/fwm.py` and `tests/test_fwm_measured_gain_integration.py`. Source/tests were frozen before the parent launched full regression. Full-suite and Grand Challenge results are recorded by the parent separately.
