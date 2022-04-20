"""
This script checks 'fixed' (stored) values of dG fitting and compares to new results generated in
`generate_test_fit_results.py'
"""

import pandas as pd
from pathlib import Path
from pyhdx.fileIO import csv_to_dataframe

test_data_dir = Path(__file__).parent / "test_data"

comparisons = {}

new_result = csv_to_dataframe(test_data_dir / "output" / "ecSecB_torch_fit.csv")
fixed_result = csv_to_dataframe(test_data_dir / "ecSecB_torch_fit_fixed.csv").rename(
    columns={"deltaG": "dG"}
)
comparisons["single_fit"] = (new_result, fixed_result)

new_result = csv_to_dataframe(
    test_data_dir / "output" / "ecSecB_torch_fit_epochs_20000.csv"
)
fixed_result = csv_to_dataframe(
    test_data_dir / "ecSecB_torch_fit_epochs_20000_fixed.csv"
).rename(columns={"deltaG": "dG"})
comparisons["single_20k"] = (new_result, fixed_result)

new_result = csv_to_dataframe(test_data_dir / "output" / "ecSecB_batch.csv")
fixed_result = csv_to_dataframe(test_data_dir / "ecSecB_batch_fixed.csv").rename(
    columns={"deltaG": "dG"}
)
comparisons["batch"] = (new_result, fixed_result)

new_result = csv_to_dataframe(test_data_dir / "output" / "ecSecB_batch_aligned.csv")
fixed_result = csv_to_dataframe(
    test_data_dir / "ecSecB_batch_aligned_fixed.csv"
).rename(columns={"deltaG": "dG"})
comparisons["batch_aligned"] = (new_result, fixed_result)


for k, (new_result, fixed_result) in comparisons.items():
    print(f"Checking {k}")
    bools = new_result.columns.get_level_values(new_result.columns.nlevels - 1) == "dG"
    new = new_result.xs("dG", level=-1, axis=1).dropna(how="any")

    bools = (
        fixed_result.columns.get_level_values(fixed_result.columns.nlevels - 1) == "dG"
    )
    old = fixed_result.xs("dG", level=-1, axis=1).dropna(how="any")

    abs_diffs = (new - old).abs()
    rel_diffs = ((new - old).abs() / old) * 100

    if new.equals(old):
        print("No differences")
    else:
        print(f"Sizes: {len(new)} ({len(old)}):")
        print(f"Differences")
        for s in ["mean", "median", "min", "max", "std"]:
            value = getattr(abs_diffs, s)()
            perc = getattr(rel_diffs, s)()
            if isinstance(value, pd.Series):
                vs = ", ".join([f"{v:.2f}" for v in value])
                ps = ", ".join([f"{p:.1f}%" for p in perc])
                print(f"{s}: {vs} ({ps})")
            else:
                print(f"{s}: {value:.2f} ({perc:.2f})%")
