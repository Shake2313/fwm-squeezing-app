import Std

set_option autoImplicit false

/-!
# GABES one-session Lean pilot (not adopted into GC)

Exact rational, scalar/isotropic covariance algebra. hbar = 1; vacuum v = 1/2.
g is the covariance/power multiplier, NOT the field-amplitude multiplier.
Independent input/reservoir noise and independent stage reservoirs are assumed.
A physical pure-loss channel has 0 <= eta <= 1. Algebraic identities below hold
for all rational eta and do not prove general real-matrix complete positivity.
-/
namespace GabesPilot

structure Channel where
  gain : Rat
  noise : Rat
  deriving Repr, DecidableEq

def run (c : Channel) (v : Rat) : Rat := c.gain * v + c.noise

/-- "after" acts on both the signal and the noise produced by "before". -/
def compose (after before : Channel) : Channel :=
  ⟨after.gain * before.gain, after.gain * before.noise + after.noise⟩

def identity : Channel := ⟨1, 0⟩
def loss (eta : Rat) : Channel := ⟨eta, (1 - eta) / 2⟩
def amplifier (g : Rat) : Channel := ⟨g, (g - 1) / 2⟩

theorem composition_action (a b : Channel) (v : Rat) :
    run (compose b a) v = run b (run a v) := by
  simp only [run, compose]
  grind

theorem composition_associative (a b c : Channel) :
    compose c (compose b a) = compose (compose c b) a := by
  cases a; cases b; cases c
  simp only [compose, Channel.mk.injEq]
  constructor <;> grind

theorem identity_after (a : Channel) : compose identity a = a := by
  cases a
  simp [compose, identity] <;> grind

theorem identity_before (a : Channel) : compose a identity = a := by
  cases a
  simp [compose, identity] <;> grind

theorem loss_preserves_vacuum (eta : Rat) :
    run (loss eta) (1 / 2) = 1 / 2 := by
  simp only [run, loss]
  grind

theorem loss_composition (eta1 eta2 : Rat) :
    compose (loss eta2) (loss eta1) = loss (eta2 * eta1) := by
  simp only [compose, loss, Channel.mk.injEq]
  constructor <;> grind

/-- Exact physical fixture: two independent 50% power-transmission losses. -/
theorem two_half_losses :
    compose (loss (1 / 2)) (loss (1 / 2)) = ⟨1 / 4, 3 / 8⟩ := by decide +kernel

theorem two_half_losses_vacuum :
    run (compose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) = 1 / 2 := by decide +kernel

/-- Incorrect: adds source noises without transmitting the earlier noise. -/
def naiveCompose (after before : Channel) : Channel :=
  ⟨after.gain * before.gain, before.noise + after.noise⟩

/-- Incorrect: drops the first reservoir altogether. -/
def missingCompose (after before : Channel) : Channel :=
  ⟨after.gain * before.gain, after.noise⟩

theorem naive_noise_witness :
    run (naiveCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) = 5 / 8 := by decide +kernel

theorem missing_noise_witness :
    run (missingCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) = 3 / 8 := by decide +kernel

theorem naive_noise_is_wrong :
    run (naiveCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) ≠ 1 / 2 := by decide +kernel

theorem missing_noise_is_wrong :
    run (missingCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) ≠ 1 / 2 := by decide +kernel

/-- Both orders have power gain 1; their added noise differs. -/
theorem equal_gain_different_noise :
    (compose (amplifier 2) (loss (1 / 2))).gain =
      (compose (loss (1 / 2)) (amplifier 2)).gain ∧
    (compose (amplifier 2) (loss (1 / 2))).noise ≠
      (compose (loss (1 / 2)) (amplifier 2)).noise := by decide +kernel

#print axioms composition_action
#print axioms composition_associative
#print axioms identity_after
#print axioms identity_before
#print axioms loss_preserves_vacuum
#print axioms loss_composition
#print axioms two_half_losses
#print axioms two_half_losses_vacuum
#print axioms naive_noise_witness
#print axioms missing_noise_witness
#print axioms naive_noise_is_wrong
#print axioms missing_noise_is_wrong
#print axioms equal_gain_different_noise

end GabesPilot
