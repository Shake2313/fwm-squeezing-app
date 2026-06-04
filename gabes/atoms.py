"""
Atom / level-scheme models and a small registry.

`AtomModel` is pure data describing a level scheme: dimension, which levels are
ground / excited, spontaneous-emission channels, dephasing channels, and which
diagonal entries shift with velocity (Doppler). It precomputes the Lindblad
super-operator and the velocity-shift super-operator once at construction.

Today only the 85Rb D1 double-Λ 4-level model is registered (exactly the level
structure the original fwm_obe.py hard-coded). Scalar 3-level schemes and Zeeman
manifolds get added in later phases.
"""
from dataclasses import dataclass, field

import numpy as np

from . import constants, hyperfine
from .core import comm_super


@dataclass
class AtomModel:
    """A level scheme. Indices are 0-based into the n_levels Hilbert space."""
    name: str
    n_levels: int
    labels: tuple
    ground: tuple
    excited: tuple
    decay: tuple          # spontaneous emission: (from_idx, to_idx, rate) tuples
    dephasing: tuple      # coherence decay on ρ_{ij}: (i, j, rate) tuples
    doppler_levels: tuple  # diagonal entries shifted by Δ_eff (= excited states)

    lindblad: np.ndarray = field(init=False, repr=False)
    S_v: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.lindblad = self._build_lindblad()
        self.S_v = self._build_velocity_shift()

    @property
    def rho_dim(self):
        return self.n_levels * self.n_levels

    def rho_index(self, row, col):
        return row * self.n_levels + col

    def _build_lindblad(self):
        """Dissipator super-operator: spontaneous emission + ground dephasing.

        Reproduces the original `_build_lindblad_fixed` for the 4-level model.
        """
        n = self.n_levels
        eye = np.eye(n, dtype=complex)
        kets = np.eye(n, dtype=complex)
        D = np.zeros((n * n, n * n), dtype=complex)

        for (i, j, rate) in self.decay:                  # decay i -> j
            L = np.sqrt(rate) * np.outer(kets[j], kets[i])   # |j⟩⟨i|
            LdL = L.conj().T @ L
            D += (np.kron(L, L.conj())
                  - 0.5 * np.kron(LdL, eye)
                  - 0.5 * np.kron(eye, LdL.T))

        for (i, j, rate) in self.dephasing:
            idx = self.rho_index(i, j)
            D[idx, idx] -= rate
        return D

    def _build_velocity_shift(self):
        """S_v = comm_super(diag over doppler_levels), so L₀(Δ_eff)=L₀(0)−Δ_eff·S_v."""
        Dee = np.zeros((self.n_levels, self.n_levels), dtype=complex)
        for k in self.doppler_levels:
            Dee[k, k] = 1.0
        return comm_super(Dee)


def _double_lambda_rb85():
    """85Rb D1 4-level double-Λ: g₁=F2, g₂=F3, e₂=F'2, e₃=F'3 (indices 0,1,2,3)."""
    ground = (0, 1)
    excited = (2, 3)
    ground_F = {0: 2, 1: 3}
    excited_F = {2: 2, 3: 3}
    decay = []
    for e in excited:
        weights = {g: hyperfine.CF2[(ground_F[g], excited_F[e])] for g in ground}
        total = sum(weights.values())
        decay.extend((e, g, constants.GAMMA * weights[g] / total) for g in ground)
    dephasing = ((0, 1, constants.GAMMA_GG), (1, 0, constants.GAMMA_GG))
    return AtomModel(
        name="double_lambda_rb85",
        n_levels=4,
        labels=("F=2", "F=3", "F'=2", "F'=3"),
        ground=ground,
        excited=excited,
        decay=tuple(decay),
        dephasing=dephasing,
        doppler_levels=excited,
    )


def two_level(gamma=None):
    """Closed 2-level absorber: g(0) ↔ e(1). Single decay channel at full Γ."""
    gamma = constants.GAMMA if gamma is None else gamma
    return AtomModel(
        name="two_level", n_levels=2, labels=("g", "e"),
        ground=(0,), excited=(1,),
        decay=((1, 0, gamma),), dephasing=(), doppler_levels=(1,),
    )


def lambda3(gamma_gg=None, gamma=None):
    """
    3-level Λ: ground g₁(0), g₂(1), excited e(2). Excited decays to each ground
    at Γ/2 (symmetric branching, total Γ). `gamma_gg` is the ground-coherence
    decay that sets the EIT/CPT (dark-resonance) linewidth — exposed as a knob.
    """
    gamma = constants.GAMMA if gamma is None else gamma
    gamma_gg = constants.GAMMA_GG if gamma_gg is None else gamma_gg
    return AtomModel(
        name="lambda3", n_levels=3, labels=("g1", "g2", "e"),
        ground=(0, 1), excited=(2,),
        decay=((2, 0, gamma / 2), (2, 1, gamma / 2)),
        dephasing=((0, 1, gamma_gg), (1, 0, gamma_gg)),
        doppler_levels=(2,),
    )


def sas_atom(n_excited=1, gamma=None):
    """
    Saturated-absorption manifold: one ground g(0) and `n_excited` excited states
    (1, 2, …) sharing it, each decaying to g at full Γ. Excited splittings are
    applied in the Hamiltonian by the scheme; all excited levels are Doppler-shifted.
    """
    gamma = constants.GAMMA if gamma is None else gamma
    n = 1 + n_excited
    excited = tuple(range(1, n))
    labels = ("g",) + tuple(f"e{i}" for i in range(1, n))
    decay = tuple((e, 0, gamma) for e in excited)
    return AtomModel(
        name=f"sas_{n_excited}", n_levels=n, labels=labels,
        ground=(0,), excited=excited, decay=decay, dephasing=(),
        doppler_levels=excited,
    )


REGISTRY = {
    "double_lambda_rb85": _double_lambda_rb85(),
}


def get(name):
    return REGISTRY[name]


def rb85_density(T):
    """
    Rb-85 atomic number density [/m³] from Steck Rb vapor pressure
    (liquid phase, T > 312.46 K). 85-Rb isotopic abundance applied.
    """
    log10_P_torr = (15.88253 - 4529.635 / T
                    + 0.00058663 * T - 2.99138 * np.log10(T))
    P_pa = 10 ** log10_P_torr * 133.322387415
    N_total = P_pa / (constants.KB * T)
    return constants.RB85_ABUNDANCE * N_total
