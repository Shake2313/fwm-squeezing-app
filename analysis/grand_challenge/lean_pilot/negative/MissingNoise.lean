import Std
import Main
open GabesPilot

-- EXPECTED REJECTION: the left side is 3/8.
theorem false_missing_vacuum :
    run (missingCompose (loss (1 / 2)) (loss (1 / 2))) (1 / 2) = 1 / 2 := by
  rw [missing_noise_witness]
  grind
