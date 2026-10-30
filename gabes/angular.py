"""
Angular-momentum algebra shared by every level-scheme builder.

Hand-rolled Racah formulas (no sympy dependency) for Clebsch-Gordan
coefficients and Wigner 3j/6j symbols; integer and half-integer arguments.
This is a leaf module (stdlib only) so atomic data (`species`), its ⁸⁵Rb D1
view (`hyperfine`) and the Zeeman builders (`zeeman`) can all import it
without an import cycle.
"""
import math


def clebsch_gordan(j1, m1, j2, m2, J, M):
    """⟨j1 m1 j2 m2 | J M⟩ via the Racah formula. Works for integer/half-integer."""
    if m1 + m2 != M:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0

    f = math.factorial
    pref = math.sqrt(
        (2 * J + 1)
        * f(int(j1 + j2 - J)) * f(int(j1 - j2 + J)) * f(int(-j1 + j2 + J))
        / f(int(j1 + j2 + J + 1))
    )
    pref *= math.sqrt(
        f(int(J + M)) * f(int(J - M))
        * f(int(j1 - m1)) * f(int(j1 + m1))
        * f(int(j2 - m2)) * f(int(j2 + m2))
    )

    kmin = int(max(0, j2 - m1 - J, j1 + m2 - J))
    kmax = int(min(j1 + j2 - J, j1 - m1, j2 + m2))
    s = 0.0
    for k in range(kmin, kmax + 1):
        denom = (f(k) * f(int(j1 + j2 - J - k)) * f(int(j1 - m1 - k))
                 * f(int(j2 + m2 - k)) * f(int(J - j2 + m1 + k)) * f(int(J - j1 - m2 + k)))
        s += (-1) ** k / denom
    return pref * s


def wigner3j(j1, j2, j3, m1, m2, m3):
    """(j1 j2 j3 ; m1 m2 m3) from the Clebsch-Gordan coefficient. 0 if Σm ≠ 0."""
    if abs(m1 + m2 + m3) > 1e-9:
        return 0.0
    return ((-1) ** int(round(j1 - j2 - m3)) / math.sqrt(2 * j3 + 1)
            * clebsch_gordan(j1, m1, j2, m2, j3, -m3))


def _fac(n):
    return math.factorial(int(round(n)))


def _is_int(x):
    return abs(x - round(x)) < 1e-9


def _tri_ok(a, b, c):
    """Triangle condition |a−b| ≤ c ≤ a+b with integer perimeter."""
    return (c >= abs(a - b) - 1e-9 and c <= a + b + 1e-9
            and _is_int(a + b + c))


def _delta(a, b, c):
    return math.sqrt(_fac(a + b - c) * _fac(a - b + c) * _fac(-a + b + c)
                     / _fac(a + b + c + 1))


def wigner6j(j1, j2, j3, j4, j5, j6):
    """{j1 j2 j3 ; j4 j5 j6} via the Racah W-coefficient sum. 0 if non-physical."""
    for (a, b, c) in ((j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3)):
        if not _tri_ok(a, b, c):
            return 0.0
    pref = _delta(j1, j2, j3) * _delta(j1, j5, j6) * _delta(j4, j2, j6) * _delta(j4, j5, j3)
    t_lo = max(j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3)
    t_hi = min(j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4)
    s = 0.0
    t = int(math.ceil(t_lo - 1e-9))
    while t <= int(math.floor(t_hi + 1e-9)):
        denom = (_fac(t - j1 - j2 - j3) * _fac(t - j1 - j5 - j6)
                 * _fac(t - j4 - j2 - j6) * _fac(t - j4 - j5 - j3)
                 * _fac(j1 + j2 + j4 + j5 - t) * _fac(j2 + j3 + j5 + j6 - t)
                 * _fac(j3 + j1 + j6 + j4 - t))
        s += (-1) ** t * _fac(t + 1) / denom
        t += 1
    return pref * s
