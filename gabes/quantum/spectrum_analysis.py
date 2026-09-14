"""Connected RF squeezing bands with explicit data gaps and scan censoring.

This layer consumes a declared linear PSD/SQL ratio. It does not construct a
source covariance, infer bandwidth from gain, or convolve an analyser filter.
Crossing precision and validity of the supplied covariance are separate facts.
"""

from dataclasses import dataclass

import numpy as np

from .contracts import AnalysisFrequencyAxis
from .readout import IntensityDifferenceSpectrum


def _owned_vector(value, size, *, boolean=False):
    raw = np.asarray(value)
    if raw.shape != (size,):
        raise ValueError("one value per RF sample required")
    if boolean:
        if raw.dtype.kind != "b":
            raise ValueError("validity and convergence masks must be boolean")
        result = np.array(raw, dtype=bool, copy=True)
    else:
        if np.iscomplexobj(raw) and np.any(raw.imag != 0):
            raise ValueError("real PSD or ratio values required")
        result = np.array(raw.real, dtype=float, copy=True)
    return np.frombuffer(result.tobytes(), dtype=result.dtype)


@dataclass(frozen=True)
class SpectralTrace:
    """Pointwise positive-RF linear ratios, with independently declared status.

    NaN/Inf samples and ``valid=False`` samples create unavailable gaps.
    ``converged=False`` separately identifies rejected numerical samples. No
    interpolation crosses either kind of gap. ``layer`` locates the supplied
    ratio at the source or detector; it never promotes its physical provenance.
    """

    analysis_axis: AnalysisFrequencyAxis
    ratio: np.ndarray
    layer: str
    provenance: str
    valid: np.ndarray | None = None
    converged: np.ndarray | None = None
    measurement_filter_scope: str = "pointwise PSD; no RBW/VBW convolution"

    def __post_init__(self):
        if not isinstance(self.analysis_axis, AnalysisFrequencyAxis):
            raise TypeError("a laboratory AnalysisFrequencyAxis is required")
        if np.any(self.analysis_axis.omega_rad_s <= 0):
            raise ValueError("strictly positive laboratory RF samples required")
        if self.layer not in {"source", "detected"}:
            raise ValueError("layer must be source or detected")
        for name in ("provenance", "measurement_filter_scope"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"explicit {name} required")
        size = len(self.analysis_axis.omega_rad_s)
        ratio = _owned_vector(self.ratio, size)
        if np.any(np.isfinite(ratio) & (ratio < 0)):
            raise ValueError("PSD/SQL ratios cannot be negative; no clipping applied")
        object.__setattr__(self, "ratio", ratio)
        for name in ("valid", "converged"):
            value = getattr(self, name)
            if value is None:
                value = np.ones(size, dtype=bool)
            object.__setattr__(self, name, _owned_vector(value, size, boolean=True))

    @classmethod
    def from_psds(cls, analysis_axis, numerator_A2_Hz, sql_A2_Hz, *, layer,
                  provenance, valid=None, converged=None,
                  measurement_filter_scope="pointwise PSD; no RBW/VBW convolution"):
        """Divide like one-sided A²/Hz PSDs; a finite SQL must be positive.

        These named arguments declare identical units and measurement planes;
        this adapter cannot establish that externally supplied traces actually
        used the same currents, balance, filter or dark-noise convention.
        Nonfinite inputs remain unavailable instead of becoming zero noise.
        """
        if not isinstance(analysis_axis, AnalysisFrequencyAxis):
            raise TypeError("a laboratory AnalysisFrequencyAxis is required")
        size = len(analysis_axis.omega_rad_s)
        numerator = _owned_vector(numerator_A2_Hz, size)
        sql = _owned_vector(sql_A2_Hz, size)
        if np.any(np.isfinite(numerator) & (numerator < 0)):
            raise ValueError("physical numerator PSD must be nonnegative")
        if np.any(np.isfinite(sql) & (sql <= 0)):
            raise ValueError("finite SQL PSD must be strictly positive")
        available = np.isfinite(numerator) & np.isfinite(sql)
        ratio = np.full(size, np.nan)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            np.divide(numerator, sql, out=ratio, where=available)
        return cls(analysis_axis, ratio, layer, provenance, valid, converged,
                   measurement_filter_scope)

    @classmethod
    def from_intensity_difference(cls, spectrum, *, contribution="total",
                                  layer="detected", provenance=None,
                                  valid=None, converged=None):
        """Adapt the established readout, retaining its quantum/electronics split.

        Use ``layer='source'`` only when the supplied readout was explicitly
        evaluated at that plane (for example with an ideal detector). The
        adapter cannot undo a detector response or optical loss.
        """
        if not isinstance(spectrum, IntensityDifferenceSpectrum):
            raise TypeError("IntensityDifferenceSpectrum required")
        if contribution not in {"quantum", "total"}:
            raise ValueError("contribution must be quantum or total")
        size = len(spectrum.analysis_axis.omega_rad_s)
        quantum = _owned_vector(spectrum.quantum_psd_A2_Hz, size)
        electronics = _owned_vector(spectrum.electronics_psd_A2_Hz, size)
        if np.any(np.isfinite(quantum) & (quantum < 0)):
            raise ValueError("physical quantum PSD must be nonnegative")
        if np.any(np.isfinite(electronics) & (electronics < 0)):
            raise ValueError("independent electronics PSD cannot be negative")
        numerator = quantum if contribution == "quantum" else quantum + electronics
        return cls.from_psds(spectrum.analysis_axis, numerator,
                             spectrum.sql_psd_A2_Hz, layer=layer,
                             provenance=(spectrum.source if provenance is None else provenance)
                             + f"; {contribution} intensity-difference PSD",
                             valid=valid, converged=converged)


@dataclass(frozen=True)
class BandEdge:
    frequency_hz: float
    kind: str                          # crossing | domain | invalid | unconverged
    bracket_hz: tuple[float, float] | None
    refinement_status: str             # interpolated | exact_sample | refined | unresolved | not_applicable
    iterations: int = 0

    @property
    def bracket_width_hz(self):
        return None if self.bracket_hz is None else self.bracket_hz[1] - self.bracket_hz[0]


@dataclass(frozen=True)
class SpectralBand:
    lower: BandEdge
    upper: BandEdge
    minimum_ratio: float                # sampled minimum only
    minimum_frequency_hz: float
    sample_count: int

    @property
    def width_hz(self):
        """A finite crossing-to-crossing width; None for every censored band."""
        if self.lower.kind == self.upper.kind == "crossing":
            return self.upper.frequency_hz - self.lower.frequency_hz
        return None

    @property
    def observed_span_hz(self):
        return self.upper.frequency_hz - self.lower.frequency_hz


@dataclass(frozen=True)
class SpectrumBandAnalysis:
    status: str                        # bands_found | no_squeezing | no_target_band | incomplete | unavailable
    bands: tuple[SpectralBand, ...]
    threshold_db: float
    threshold_ratio: float
    minimum_ratio: float | None
    minimum_frequency_hz: float | None
    invalid_sample_count: int
    unconverged_sample_count: int
    complete_domain: bool               # all sampled values accepted, not infinite RF coverage
    layer: str
    provenance: str
    measurement_filter_scope: str
    rf_refinement_status: str           # sampled_only | refined | unresolved | no_crossings
    crossing_tolerance_hz: float | None
    scope: str = ("strict below-threshold components on the sampled RF domain; "
                  "linear-ratio interpolation; sampled minima; unresolved bands "
                  "between samples require independent RF-grid refinement")


def _evaluate_ratio(evaluator, frequency_hz):
    value = np.asarray(evaluator(float(frequency_hz)))
    if value.ndim != 0 or np.iscomplexobj(value):
        raise ValueError("crossing evaluator must return one real linear ratio")
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError("crossing evaluator returned unavailable or nonphysical data")
    return value


def _crossing_edge(f0, r0, f1, r1, threshold, evaluator, tolerance, max_iterations):
    if evaluator is not None:
        for frequency, sampled in ((f0, r0), (f1, r1)):
            evaluated = _evaluate_ratio(evaluator, frequency)
            if not np.isclose(evaluated, sampled, rtol=1e-10, atol=1e-14):
                raise ValueError("crossing evaluator disagrees with the supplied trace endpoints")
            if np.sign(evaluated-threshold) != np.sign(sampled-threshold):
                raise ValueError("crossing evaluator changes the sampled threshold bracket")
    if r0 == threshold:
        return BandEdge(f0, "crossing", (f0, f0), "exact_sample")
    if r1 == threshold:
        return BandEdge(f1, "crossing", (f1, f1), "exact_sample")
    # Form the bounded dimensionless fraction first. Multiplying a large but
    # finite PSD ratio by an RF span can overflow even for a valid crossing.
    fraction = (threshold-r0)/(r1-r0)
    estimate = f0 + fraction*(f1-f0)
    if evaluator is None:
        return BandEdge(float(estimate), "crossing", (f0, f1), "interpolated")
    lo, hi, ylo, yhi = f0, f1, r0-threshold, r1-threshold
    iterations = 0
    while hi-lo > tolerance and iterations < max_iterations:
        mid = lo + (hi-lo)/2
        if mid in (lo, hi):
            break
        value = _evaluate_ratio(evaluator, mid)
        ym = value-threshold
        iterations += 1
        if ym == 0:
            lo = hi = mid
            ylo = yhi = 0.
            break
        if (ym < 0) == (ylo < 0):
            lo, ylo = mid, ym
        else:
            hi, yhi = mid, ym
    # Use the bracket midpoint so its maximum crossing-location error is
    # explicitly bounded by half the returned bracket, subject to continuity.
    status = "refined" if hi-lo <= tolerance else "unresolved"
    return BandEdge(float(lo+(hi-lo)/2), "crossing", (float(lo), float(hi)), status, iterations)


def analyze_squeezing_bands(trace, *, threshold_db=0., crossing_evaluator=None,
                            crossing_tolerance_hz=None, max_refinement_iterations=60):
    """Extract every strict R(f)<10**(threshold_db/10) connected component.

    A finite width requires two observed crossings. Scan edges and invalid or
    unconverged holes remain censored. Equality samples/plateaus separate strict
    components, even when below-threshold samples lie on both sides.

    Optional local bisection evaluates the *same* continuous provider in each
    observed crossing bracket. ``crossing_evaluator(frequency_hz)`` returns its
    validated linear ratio; unavailable evaluations fail explicitly. The
    provider must match both sampled endpoints (rtol=1e-10, atol=1e-14) without
    changing either endpoint's strict threshold side or equality. Refinement
    does not certify global RF coverage, covariance convergence or
    model/measurement uncertainty and cannot discover entirely missed bands.
    """
    if not isinstance(trace, SpectralTrace):
        raise TypeError("SpectralTrace required")
    threshold_db = float(threshold_db)
    if not np.isfinite(threshold_db) or threshold_db > 0:
        raise ValueError("finite SQL or below-SQL target threshold in dB required")
    threshold = float(10.**(threshold_db/10))
    if threshold == 0:
        raise ValueError("target threshold underflows the linear-ratio representation")
    if crossing_evaluator is not None:
        if not callable(crossing_evaluator):
            raise TypeError("crossing evaluator must be callable")
        if crossing_tolerance_hz is None:
            raise ValueError("declare a crossing tolerance in Hz for local refinement")
        crossing_tolerance_hz = float(crossing_tolerance_hz)
        if not np.isfinite(crossing_tolerance_hz) or crossing_tolerance_hz <= 0:
            raise ValueError("positive finite crossing tolerance required")
    elif crossing_tolerance_hz is not None:
        raise ValueError("crossing tolerance requires a provider for refinement")
    if (isinstance(max_refinement_iterations, bool)
            or not isinstance(max_refinement_iterations, (int, np.integer))
            or max_refinement_iterations < 1):
        raise ValueError("positive integer maximum refinement iterations required")
    frequency, ratio = trace.analysis_axis.frequency_hz, trace.ratio
    valid = trace.valid & np.isfinite(ratio)
    usable = valid & trace.converged
    below = usable & (ratio < threshold)

    def edge(index, step):
        neighbor = index+step
        if not 0 <= neighbor < len(ratio):
            return BandEdge(float(frequency[index]), "domain", None, "not_applicable")
        if not usable[neighbor]:
            kind = "invalid" if not valid[neighbor] else "unconverged"
            return BandEdge(float(frequency[index]), kind, None, "not_applicable")
        left, right = sorted((index, neighbor))
        return _crossing_edge(float(frequency[left]), float(ratio[left]),
                              float(frequency[right]), float(ratio[right]),
                              threshold, crossing_evaluator, crossing_tolerance_hz,
                              max_refinement_iterations)

    bands = []
    starts = np.flatnonzero(below & ~np.r_[False, below[:-1]])
    ends = np.flatnonzero(below & ~np.r_[below[1:], False])
    for start, end in zip(starts, ends):
        minimum = int(start + np.argmin(ratio[start:end+1]))
        bands.append(SpectralBand(edge(int(start), -1), edge(int(end), 1),
                                  float(ratio[minimum]), float(frequency[minimum]),
                                  int(end-start+1)))
    complete = bool(np.all(usable))
    if not np.any(usable):
        status, minimum_ratio, minimum_frequency = "unavailable", None, None
    else:
        minimum = int(np.flatnonzero(usable)[np.argmin(ratio[usable])])
        minimum_ratio, minimum_frequency = float(ratio[minimum]), float(frequency[minimum])
        absence = "no_squeezing" if threshold_db == 0 else "no_target_band"
        status = "bands_found" if bands else (absence if complete else "incomplete")
    crossing_statuses = [edge.refinement_status for band in bands
                         for edge in (band.lower, band.upper) if edge.kind == "crossing"]
    if not crossing_statuses:
        refinement = "no_crossings"
    elif "unresolved" in crossing_statuses:
        refinement = "unresolved"
    else:
        refinement = "sampled_only" if crossing_evaluator is None else "refined"
    return SpectrumBandAnalysis(status, tuple(bands), threshold_db, threshold,
                                minimum_ratio, minimum_frequency,
                                int(np.count_nonzero(~valid)),
                                int(np.count_nonzero(valid & ~trace.converged)),
                                complete, trace.layer, trace.provenance,
                                trace.measurement_filter_scope, refinement,
                                crossing_tolerance_hz)
