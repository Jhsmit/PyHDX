"""Load HDX-MS data from yaml spec and perform initial guess of exchange rates"""
from pyhdx.batch_processing import load_from_yaml
from pathlib import Path
from pyhdx.fitting import fit_rates_weighted_average, fit_rates_half_time_interpolate
import yaml
from pyhdx.local_cluster import default_client

current_dir = Path(__file__).parent
output_dir = current_dir / 'guesses'
output_dir.mkdir(exist_ok=True)
data_dir = current_dir.parent / 'tests' / 'test_data'
yaml_stream = Path(current_dir / 'yaml_files' / 'SecB.yaml').read_text()
data_dict = yaml.safe_load(yaml_stream)

# Requires local_cluster.py to be running (or other Dask client on default address in config)
client = default_client()


for name, dic in data_dict.items():
    print(name)
    dic = data_dict[name]
    hdxm = load_from_yaml(dic, data_dir=data_dir)

    #Uncomment/adjust path to save sequence info + intrinsic rates
    #hdxm.coverage.protein.to_file(f'{name}_sequence_info.txt', fmt='pprint')

    fr = fit_rates_weighted_average(hdxm, client=client)
    fr.output.to_file(output_dir / f'{name}_rates_guess.csv')



