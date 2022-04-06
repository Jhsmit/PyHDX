"""
Execute a PyHDX data processing pipeline according to a yaml jobfile specification


"""


from pathlib import Path
from pyhdx.batch_processing import JobParser
import yaml

#%%
# Pycharm scientific mode
if '__file__' not in locals():
    __file__ = Path().cwd() / 'templates' / 'script.py'

current_dir = Path(__file__).parent
output_dir = current_dir / 'output'
output_dir.mkdir(exist_ok=True)
test_data_dir = current_dir.parent / 'tests' / 'test_data'
input_dir = test_data_dir / 'input'

#%%

job_spec = yaml.safe_load((input_dir / 'jobfile.yaml').read_text())
job_parser = JobParser(job_spec, cwd=input_dir)
job_parser.execute()