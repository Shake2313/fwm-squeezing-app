"""Independent one-period density-operator QRT, with an exact geometric tail.

No atomic drift, diffusion, Floquet mean coefficients or lifted solver imports.
Integrates the full time-dependent Liouvillian, finds its periodic fixed point,
and resums all subsequent periods for each positive-delay regression integral.
"""

import numpy as np
from scipy.integrate import solve_ivp


def time_domain_qrt(l0, lplus, lminus, beat_rad_s, operators, base_frequencies_rad_s,
                    output_harmonics, *, phase_samples=32, rtol=2e-10, atol=2e-12):
    l0, lp, lm, ops = [np.asarray(x, complex) for x in (l0, lplus, lminus, operators)]
    beat = float(beat_rad_s)
    frequencies = np.asarray(base_frequencies_rad_s, float)
    harmonics = np.asarray(output_harmonics, int)
    n = ops.shape[-1]
    size, m = n*n, len(ops)
    if (l0.shape != (size, size) or lp.shape != l0.shape or lm.shape != l0.shape
            or ops.shape != (m, n, n) or not np.isfinite(beat) or beat == 0
            or phase_samples < 2*max(abs(harmonics))+3):
        raise ValueError("matched generators, operators, beat and sufficient reference phases required")
    speed, sign = abs(beat), np.sign(beat)
    def generator(theta):
        return l0+lp*np.exp(-1j*sign*theta)+lm*np.exp(1j*sign*theta)
    eye = np.eye(size, dtype=complex)
    flow = solve_ivp(lambda theta, y: (generator(theta)@y.reshape(size, size)/speed).reshape(-1),
        (0., 2*np.pi), eye.reshape(-1), method="DOP853", dense_output=True, rtol=rtol, atol=atol)
    if not flow.success:
        raise ValueError("reference period propagation failed")
    period_map = flow.y[:, -1].reshape(size, size)
    trace = np.eye(n).reshape(-1)
    fixed = period_map-eye
    fixed[0] = trace
    rhs = np.zeros(size, complex); rhs[0] = 1
    rho0 = np.linalg.solve(fixed, rhs)
    phase_grid = 2*np.pi*np.arange(phase_samples)/phase_samples
    rho = np.array([(flow.sol(theta).reshape(size, size)@rho0).reshape(n, n) for theta in phase_grid])
    actual_omega = (frequencies[:, None]+harmonics[None, :]*beat).reshape(-1)
    nh, nf = len(harmonics), len(frequencies)
    positive = np.zeros((nf, nh, m, nh, m), complex)
    readout = ops.swapaxes(-1, -2).reshape(m, size)
    maximum_regression_trace = 0.
    for theta0, state in zip(phase_grid, rho):
        # Dimensionless one-period integral Q; multiply by 1/abs(beat) below.
        initial = np.zeros((1+len(actual_omega), size, size), complex)
        initial[0] = eye
        def ode(theta, raw):
            values = raw.reshape(initial.shape)
            result = np.empty_like(values)
            result[0] = generator(theta0+theta)@values[0]/speed
            result[1:] = np.exp(1j*actual_omega*theta/speed)[:, None, None]*values[0]
            return result.reshape(-1)
        solution = solve_ivp(ode, (0., 2*np.pi), initial.reshape(-1), method="DOP853", rtol=rtol, atol=atol)
        if not solution.success:
            raise ValueError("reference regression period integration failed")
        values = solution.y[:, -1].reshape(initial.shape)
        p, integrals = values[0], values[1:]/speed
        means = np.einsum("iab,ba->i", ops, state)
        sources = ((ops-means[:, None, None]*np.eye(n))@state).reshape(m, size).T
        for k, w in enumerate(actual_omega):
            # Border the trace mode, including exact RF=0 / integer-harmonic w.
            bordered = np.zeros((size+1, size+1), complex)
            bordered[:size, :size] = eye-np.exp(1j*w*2*np.pi/speed)*p
            bordered[:size, size], bordered[size, :size] = state.reshape(-1), trace
            source = np.vstack([sources, np.zeros((1, m))])
            tail = np.linalg.solve(bordered, source)[:size]
            maximum_regression_trace = max(maximum_regression_trace, float(np.max(abs(trace@tail))))
            integral = readout@integrals[k]@tail
            fi, hi = divmod(k, nh)
            for ki, harmonic in enumerate(harmonics):
                positive[fi, hi, :, ki, :] += integral*np.exp(1j*(harmonics[hi]-harmonic)*sign*theta0)/phase_samples
    positive = positive.reshape(nf, nh*m, nh*m)
    spectrum = positive+positive.conj().swapaxes(-1, -2)
    return {"ordered_spectrum": spectrum, "phase_states": rho,
            "physical_phases_rad": sign*phase_grid, "period_map": period_map,
            "maximum_regression_trace": maximum_regression_trace,
            "periodic_fixed_point_residual": float(np.linalg.norm((period_map-eye)@rho0)),
            "minimum_phase_state_eigenvalue": float(np.linalg.eigvalsh((rho+rho.conj().swapaxes(-1, -2))/2).min()),
            "phase_samples": phase_samples, "rtol": rtol, "atol": atol}
