"""Read-only independent ordered-QRT pilot; not imported by production.

This uses covariance initial sources and block-tridiagonal regression rather
than R D R-dagger. All source harmonics, including mean-product convolution,
are retained within the declared finite response ladder.
"""
from pathlib import Path
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gabes import constants as c, core
from gabes.schemes import fwm
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.diffusion import stationary_atomic_noise, traceless_hermitian_basis
from gabes.quantum.periodic import periodic_atomic_noise
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs, thermal_reset_channels
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.quantum.contracts import OpticalDetunings


def explicit_same_production_reservoirs(temperature_K):
    atom = fwm.collisional_atom(temperature_K)
    gg = atom.ground_collision_dephasing_rate
    optical = max((rate for i, j, rate in atom.dephasing
                   if i in fwm.GROUND_STATES and j in fwm.EXCITED_STATES), default=0.)
    if optical < gg/4:
        raise ValueError("specified coherence damping has no declared positive diagonal-jump decomposition")
    radiative = ExplicitReservoirs.from_atom(
        fwm.atoms.double_lambda_rb85(gamma_gg=0.), source="same production radiative channels")
    extra = (
        CollapseChannel("coll:ground", np.sqrt(gg/2)*np.diag([1., -1., 0., 0.]),
                        "dephasing", "exact target coherence rates"),
        CollapseChannel("coll:optical", np.sqrt(2*(optical-gg/4))*np.diag([0., 0., 1., 1.]),
                        "dephasing", "exact target coherence rates"),
    ) + thermal_reset_channels(atom.transit_reset_rate, atom.thermal_reset_populations,
                               source="same production thermal reset")
    reservoirs = ExplicitReservoirs(4, radiative.channels+extra)
    parity = float(np.linalg.norm(reservoirs.dissipator()-atom.lindblad)
                   / np.linalg.norm(atom.lindblad))
    return atom, reservoirs, parity


def block_tridiagonal_action(diagonal, lower, upper, rhs):
    """Solve block (lower, diagonal, upper) with one or many RHS columns."""
    blocks = len(diagonal)
    solved_rhs, solved_upper = [], []
    for index in range(blocks):
        center = diagonal[index].copy()
        target = rhs[index].copy()
        if index:
            center -= lower @ solved_upper[-1]
            target -= lower @ solved_rhs[-1]
        joined = np.concatenate((target, upper), axis=-1)
        result = np.linalg.solve(center, joined)
        solved_rhs.append(result[..., :target.shape[-1]])
        solved_upper.append(result[..., target.shape[-1]:])
    result = [None]*blocks
    result[-1] = solved_rhs[-1]
    for index in range(blocks-2, -1, -1):
        result[index] = solved_rhs[index]-solved_upper[index]@result[index+1]
    return np.stack(result)


def direct_periodic_qrt(h0, harmonic_drive, reservoirs, beat, state_harmonics,
                        lowering_operators, physical_frequencies, response_order):
    """Positive ordered dipole spectra via exact connected QRT source action.

    States use sum rho_n exp(-i n beat t). Frequencies are physical atomic
    frame frequencies; no quasifrequency folding is needed. A finite response
    ladder must be converged independently. Returned signed real rates are
    never clipped. This computes atomic spectra [seconds], not field rates.
    """
    start = time.perf_counter()
    n = reservoirs.n_levels
    operators = traceless_hermitian_basis(n)
    basis = operators.reshape(n*n-1, n*n).T
    states = np.asarray(state_harmonics)
    mean_order = len(states)//2
    harmonic_indices = np.arange(-mean_order, mean_order+1)
    generators = [core.comm_super(harmonic_drive.conj().T), reservoirs.generator(h0),
                  core.comm_super(harmonic_drive)]
    drift_minus, drift_zero, drift_plus = [basis.conj().T@g@basis for g in generators]
    means = np.einsum("iab,qba->qi", operators, states)
    products = operators[:, None] @ operators[None, :]
    moments = np.einsum("ijab,qba->qij", products, states)
    covariance = {}
    for h in range(-2*mean_order, 2*mean_order+1):
        item = np.zeros((n*n-1, n*n-1), complex)
        if abs(h) <= mean_order:
            item += moments[h+mean_order]
        for k in harmonic_indices:
            if abs(h-k) <= mean_order:
                item -= np.outer(means[k+mean_order], means[h-k+mean_order])
        covariance[h] = item
    ladder = np.arange(-response_order, response_order+1)
    zero = np.zeros_like(drift_zero)
    outputs = []
    for lowering, frequency in zip(lowering_operators, physical_frequencies):
        w = np.einsum("ab,iba->i", lowering, operators)
        source = np.array([np.column_stack((covariance.get(int(h), zero)@w.conj(),
                                           covariance.get(int(h), zero).T@w.conj()))
                           for h in ladder])
        diagonal = drift_zero[None] + 1j*(frequency+ladder*beat)[:, None, None]*np.eye(n*n-1)
        result = block_tridiagonal_action(diagonal, drift_plus, drift_minus, -source)
        rates = 2*np.real(w@result[response_order])
        # Two orderings differ by the dipole commutator. chi per Rabi uses /2.
        chi = -.5j*(w@(result[response_order, :, 0]-result[response_order, :, 1]))
        outputs.append({"absorption_seconds": float(rates[0]),
                        "emission_seconds": float(rates[1]),
                        "chi_per_rabi_real_s": float(chi.real),
                        "chi_per_rabi_imag_s": float(chi.imag)})
    return {"arms": outputs, "seconds": time.perf_counter()-start,
            "mean_order": mean_order, "response_order": response_order,
            "source_covariance_harmonics": 2*mean_order}


def run_pilot():
    temperature = 386.15
    atom, reservoirs, parity = explicit_same_production_reservoirs(temperature)
    op = c.rabi_freq(.2948, 530e-6)
    os = c.rabi_freq(3.918e-6, 330e-6)
    delta = 2*np.pi*23.7e6
    beat = -c.OMEGA_HF+delta
    h0 = fwm.pump_hamiltonian_at_Deff_zero(op, op)-2*np.pi*1e9*np.diag([0., 0., 1., 1.])
    lowering = []
    for ground in [fwm.G2, fwm.G1]:
        d = np.zeros((4, 4), complex)
        for e in fwm.EXCITED_STATES:
            d[ground, e] = fwm.TRANSITION_DIPOLE_SCALE[ground, e]
        lowering.append(d)
    rows = []
    for scale in [0., 1e-6, 1.]:
        v = os*np.sqrt(scale)/2*lowering[0].conj().T
        periodic = periodic_atomic_noise(h0, v, reservoirs, beat, mean_order=3, phase_samples=32)
        spectrum = periodic.lift(3).spectrum(GeneratorFrequencyAxis([0.], "pilot Floquet quasifrequency"))
        baseline = []
        for d, h in zip(lowering, [1, -1]):
            w = np.einsum("ab,iba->i", d, periodic.operators)
            baseline.append({ordering: float(np.real(w@spectrum.block(h, h, ordering=ordering)[0]@w.conj()))
                             for ordering in ["greater", "lesser"]})
        independent = [direct_periodic_qrt(h0, v, reservoirs, beat, periodic.state_harmonics,
                                          lowering, [beat, -beat], order) for order in [2, 3, 4]]
        rows.append({"seed_power_ratio": scale,
                     "mean_population": np.diag(periodic.state_harmonics[3]).real.tolist(),
                     "RDRdagger_reference": baseline, "direct_connected_QRT": independent})
    # Independent two-level Lorentzian analytic controls (not a diffusion identity).
    analytic_errors = []
    for up in [0., .5e6, 2e6]:
        down, center = 1e6, 1e7
        ds = np.array([[0., 1.], [0., 0.]])
        rs = ExplicitReservoirs(2, (
            CollapseChannel("down", np.sqrt(down)*ds, "emission", "analytic control"),
            CollapseChannel("up", np.sqrt(up)*ds.T, "incoherent_pump", "analytic control")))
        hs = np.diag([0., center])
        rho = np.diag([down/(down+up), up/(down+up)])
        harmonic = np.stack((np.zeros((2, 2)), rho, np.zeros((2, 2))))
        for frequency in [0., center, center+3e6]:
            result = direct_periodic_qrt(hs, np.zeros((2, 2)), rs, 2e7,
                                         harmonic, [ds], [frequency], 1)["arms"][0]
            width = (down+up)/2
            line = 2*width/(width**2+(frequency-center)**2)
            expected = np.array([rho[0, 0]*line, rho[1, 1]*line])
            actual = np.array([result["absorption_seconds"], result["emission_seconds"]])
            analytic_errors.append(float(np.max(np.abs(actual-expected))/np.max(expected)))
    # Same finite state as the production seed-frame solve, not an averaged rho0.
    v = os/2*lowering[0].conj().T
    states = core.floquet_solve_truncated(reservoirs.generator(h0), core.comm_super(v),
        core.comm_super(v.conj().T), beat, [0.], np.zeros((16, 16)), 4,
        n_f=4, return_harmonics=True)[0]
    seeded_h0 = fwm.static_hamiltonian_at_Deff_zero(op, op, os, delta, -1)
    cp, cm = fwm.sideband_template(op, op, 0., -1)
    seeded = core.floquet_solve_truncated(core.build_liouvillian(seeded_h0, atom), cp, cm,
        -beat, [2*np.pi*1e9], atom.S_v, 4, n_f=4, return_harmonics=True)[0]
    phases = 2*np.pi*np.arange(64)/64
    pump_states = np.einsum("tq,qij->tij", np.exp(-1j*np.outer(phases, np.arange(-4, 5))), states)
    seed_states = np.einsum("tq,qij->tij", np.exp(+1j*np.outer(phases, np.arange(-4, 5))), seeded)
    gauge = np.tile(np.eye(4, dtype=complex), (64, 1, 1))
    gauge[:, fwm.G2, fwm.G2] = np.exp(-1j*phases)
    mapped = gauge@pump_states@gauge.conj().transpose(0, 2, 1)
    gauge_error = float(np.max(np.linalg.norm(mapped-seed_states, axis=(1, 2))))
    velocity_controls = []
    for velocity in [-200., 0., 200.]:
        hv = h0+fwm.K_VEC*velocity*np.diag([0., 0., 1., 1.])
        orders = []
        for order in [3, 4]:
            state = core.floquet_solve_truncated(reservoirs.generator(hv),
                core.comm_super(v), core.comm_super(v.conj().T), beat, [0.],
                np.zeros((16, 16)), 4, n_f=order, return_harmonics=True)[0]
            orders.append(direct_periodic_qrt(hv, v, reservoirs, beat, state,
                                               lowering, [beat, -beat], order))
        rates = np.array([[arm[key] for arm in result["arms"]
                          for key in ["absorption_seconds", "emission_seconds"]]
                         for result in orders])
        velocity_controls.append({"velocity_m_s": velocity, "results": orders,
            "mean_response_order_3_vs_4_relative_error": float(np.max(
                abs(rates[1]-rates[0])/np.maximum(abs(rates[1]), np.finfo(float).tiny)))})
    result = {"origin": "AI independent ordered QRT derivation/pilot", "marker": "⟦AI⟧",
            "date": "2026-10-16", "new_atomic_velocity_classes": 3,
            "settings": {"T_K": temperature, "pump_W": .2948, "pump_waist_m": 530e-6,
                         "seed_W": 3.918e-6, "seed_waist_m": 330e-6,
                         "one_photon_detuning_GHz": 1., "internal_delta_MHz": 23.7,
                         "conjugate_carrier_input": "vacuum; pump+probe prescribed local state"},
            "explicit_dissipator_relative_parity": parity, "rows": rows,
            "two_level_analytic_cases": 9,
            "two_level_analytic_maximum_relative_error": max(analytic_errors),
            "production_finite_state_gauge_phase_count": 64,
            "production_finite_state_gauge_maximum_norm_error": gauge_error,
            "velocity_order_controls": velocity_controls,
            "scope": "read-only local atomic conditional spectra; no ensemble propagation or fit"}
    finite = fwm.chi_matrix_table(op, op, os, os, np.array([delta]),
        np.array([2*np.pi*1e9]), -1, atom=atom, n_f=3)
    result["finite_secant_control"] = {
        "chi_ss_real": float(finite[0][0, 0].real),
        "chi_ss_imag": float(finite[0][0, 0].imag),
        "chi_cc_real": float(finite[3][0, 0].real),
        "chi_cc_imag": float(finite[3][0, 0].imag),
        "scope": "finite independent probe-only/conjugate-only extraction; not one common-state Jacobian"}
    frequencies = optical_carriers(OpticalDetunings(2*np.pi*1e9, delta))[1:]
    prefactor = (fwm.hyperfine.number_density(temperature)*frequencies*c.DIPOLE_D1**2
                 *(.74/12)/(2*c.EPS_0*c.HBAR*c.C_LIGHT))
    result["field_normalization"] = {
        "prefactor": "n*omega*d^2*line_strength/(24*epsilon0*hbar*c)",
        "line_strength": .74, "macroscopic_average": 1/12,
        "origin": "inherited production weak-field normalization; not a new absolute reduction certification",
        "prefactor_per_m_per_second": prefactor.tolist()}
    for row in rows:
        rates = []
        for k, a in enumerate(row["direct_connected_QRT"][1]["arms"]):
            absorption = prefactor[k]*a["absorption_seconds"]
            emission = prefactor[k]*a["emission_seconds"]
            rates.append({"alpha_abs_m_inverse": float(absorption),
                          "alpha_em_m_inverse": float(emission),
                          "alpha_net_m_inverse": float(absorption-emission),
                          "absorption_channel_T_12p5mm": float(np.exp(-absorption*.0125))})
        row["conditional_legacy_weak_normalization_rates"] = rates
    return result


if __name__ == "__main__":
    result = run_pilot()
    target = Path(__file__).with_suffix(".json")
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=True, indent=2))
