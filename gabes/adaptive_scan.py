"""
Adaptive sampling of a fixed display scan (solve few points, interpolate the rest).

A plotted scan often has a known display grid (e.g. 401 probe detunings) while
the expensive physics is smooth almost everywhere and sharp only near a few
resonances.  :func:`refine_scan_nodes` starts from a coarse subset, predicts each
interval midpoint from the current nodes, solves the midpoint exactly, and
splits only the intervals whose prediction disagrees.  Every solved midpoint is
kept as a node, so no solve is wasted.

Interpolate a quantity that is linear in the physics (a susceptibility), not an
exponentiated observable (a gain), and judge the disagreement in the observable
the user reads.  :func:`cubic_spline` is a dependency-free not-a-knot spline with
the same end conditions as ``scipy.interpolate.CubicSpline``.
"""
import numpy as np


def cubic_spline(x_nodes, y_nodes, x):
    """Not-a-knot cubic spline of ``y_nodes`` (n, ...) evaluated at ``x``.

    Real or complex values; trailing dimensions are interpolated together.
    Fewer than four nodes fall back to piecewise-linear interpolation.
    """
    xn = np.asarray(x_nodes, dtype=float)
    yn = np.asarray(y_nodes)
    x = np.asarray(x, dtype=float)
    n = xn.size
    if n != yn.shape[0]:
        raise ValueError("x_nodes and y_nodes must have the same length")
    if n < 2 or np.any(np.diff(xn) <= 0.0):
        raise ValueError("x_nodes must be strictly increasing with at least two nodes")
    tail = yn.shape[1:]
    Y = yn.reshape(n, -1)
    if n < 4:
        cols = [np.interp(x, xn, Y[:, j].real) + (1j * np.interp(x, xn, Y[:, j].imag)
                if np.iscomplexobj(Y) else 0.0) for j in range(Y.shape[1])]
        return np.stack(cols, axis=-1).reshape(x.shape + tail)

    h = np.diff(xn)
    slope = (Y[1:] - Y[:-1]) / h[:, None]
    lower = np.zeros(n)
    diag = np.zeros(n)
    upper = np.zeros(n)
    rhs = np.zeros((n, Y.shape[1]), dtype=np.result_type(Y, float))
    d0 = xn[2] - xn[0]
    diag[0], upper[0] = h[1], d0
    rhs[0] = ((h[0] + 2.0 * d0) * h[1] * slope[0] + h[0] ** 2 * slope[1]) / d0
    lower[1:-1] = h[1:]
    diag[1:-1] = 2.0 * (h[:-1] + h[1:])
    upper[1:-1] = h[:-1]
    rhs[1:-1] = 3.0 * (h[1:, None] * slope[:-1] + h[:-1, None] * slope[1:])
    d1 = xn[-1] - xn[-3]
    lower[-1], diag[-1] = d1, h[-2]
    rhs[-1] = (h[-1] ** 2 * slope[-2] + (2.0 * d1 + h[-1]) * h[-2] * slope[-1]) / d1

    # Thomas algorithm (tridiagonal elimination); the matrix is real, the RHS may be complex.
    c = np.empty(n)
    d = np.empty_like(rhs)
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denom = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / denom
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / denom
    m = np.empty_like(rhs)
    m[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        m[i] = d[i] - c[i] * m[i + 1]

    idx = np.clip(np.searchsorted(xn, x, side="right") - 1, 0, n - 2)
    t = (x - xn[idx])[..., None]
    hi = h[idx][..., None]
    c2 = (3.0 * slope[idx] - 2.0 * m[idx] - m[idx + 1]) / hi
    c3 = (m[idx] + m[idx + 1] - 2.0 * slope[idx]) / hi ** 2
    out = Y[idx] + t * (m[idx] + t * (c2 + t * c3))          # Horner
    return out.reshape(x.shape + tail)


def refine_scan_nodes(n_points, solve, disagrees, *, start_stride=32, max_rounds=16):
    """Bisection refinement of display indices ``0..n_points-1``.

    ``solve(rows)`` computes and stores exact values for new display rows.
    ``disagrees(nodes, candidates)`` returns a boolean array: True where the
    prediction built from ``nodes`` misses the exact value at ``candidates``.
    Returns ``(sorted node indices, rounds)``; both scan ends are always nodes.
    """
    n_points = int(n_points)
    if n_points < 2:
        raise ValueError("a scan needs at least two points")
    stride = max(1, int(start_stride))
    nodes = sorted(set(range(0, n_points, stride)) | {n_points - 1})
    solve(np.asarray(nodes))
    active = [(a, b) for a, b in zip(nodes[:-1], nodes[1:]) if b - a > 1]
    rounds = 1
    while active and rounds < int(max_rounds):
        mids = [(a + b) // 2 for a, b in active]
        solve(np.asarray(mids))
        flags = np.asarray(disagrees(np.asarray(nodes), np.asarray(mids)), dtype=bool)
        nodes = sorted(set(nodes) | set(mids))
        nxt = []
        for (a, b), m, flag in zip(active, mids, flags):
            if flag:
                nxt += [iv for iv in ((a, m), (m, b)) if iv[1] - iv[0] > 1]
        active = nxt
        rounds += 1
    return np.asarray(nodes), rounds
