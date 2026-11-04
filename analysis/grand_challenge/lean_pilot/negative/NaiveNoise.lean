import Std
import Main
open GabesPilot

-- EXPECTED REJECTION: the left side is 5/8.
theorem false_naive_vacuum :
    run (naiveCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) = 1 / 2 := by
  rw [naive_noise_witness]
  grind
