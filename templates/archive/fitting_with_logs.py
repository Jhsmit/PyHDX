"""Perform fitting with a range of regularizers"""
from pyhdx.batch_processing import StateParser
from pathlib import Path
from pyhdx.fitting import fit_gibbs_global_batch
import yaml
import time
from datetime import datetime
from pyhdx.fileIO import csv_to_protein
from pyhdx import VERSION_STRING

current_dir = Path(__file__).parent
test_data_dir = current_dir.parent / "tests" / "test_data"
input_dir = test_data_dir / "input"
yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
data_dict = yaml.safe_load(yaml_stream)

output_dir = current_dir / "fit"
output_dir.mkdir(exist_ok=True)

parser = StateParser(data_dict, data_src=input_dir)
hdx_set = parser.load_hdxmset()

rates_list = [
    csv_to_protein(current_dir / "guesses" / f"{name}_rates_guess.csv")["rate"]
    for name in data_dict.keys()
]
gibbs_guess = hdx_set.guess_deltaG(rates_list)

log_file = output_dir / "fitting_log.txt"
now = datetime.now()
date = f'# {now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())})'

lines = [VERSION_STRING, date]

r2 = 0.5
for r1 in [0, 0.01, 0.25, 0.5, 1]:
    t0 = time.time()
    result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=1000, r1=r1, r2=r2)
    t1 = time.time()

    block = "--------------------------"
    regularizers = f"Regualizer 1: {r1}  Regualizer 2: {r2}"
    loss = (
        f"Total_loss {result.total_loss:.2f}, mse_loss {result.mse_loss:.2f}, reg_loss {result.reg_loss:.2f}"
        f"({result.regularization_percentage:.2f}%)"
    )
    time_elapsed = f"Time elapsed: {(t1 - t0):.2f} s"
    epochs = f"Number of epochs: {result.metadata['epochs_run']}"

    result.output.to_csv(output_dir / f"fit_output_r1_{r1}_r2_{r2}.csv")  # , na_rep='NaN')
    result.output.to_file(
        output_dir / f"fit_output_r1_{r1}_r2_{r2}.txt", fmt="pprint", na_rep="NaN"
    )

    lines += ["", block, regularizers, loss, epochs, time_elapsed, block]

log_file.write_text("\n".join(lines))
