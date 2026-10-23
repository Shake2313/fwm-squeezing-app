"""Reproduce the historical bounded reference comparison, outside the frontend."""
from pathlib import Path
import argparse
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gabes.schemes.fwm import FWMScheme


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("output/fwm_experiment_20261022/measured_reference_comparison.png"))
    args = parser.parse_args()
    scheme = FWMScheme()
    params = dict(scheme.presets()[0].values, measured_gain_calibration=True)
    comparison = scheme.internal_measured_gain_comparison()
    figure = comparison.render(comparison.compute(params), params)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=150)
    import matplotlib.pyplot as plt
    plt.close(figure)
    print(args.output)


if __name__ == "__main__":
    main()
