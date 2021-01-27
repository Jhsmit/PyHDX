from pyhdx.batch_processing import load_from_yaml, do_fitting_from_yaml
import yaml
from pathlib import Path


root_dir = Path().resolve()
data_dir = root_dir.parent / 'tests' / 'test_data'

print(data_dir)
yaml_file_path = Path('SecB.yaml')
SecB_dic = yaml.load(yaml_file_path.read_text())

kf = load_from_yaml(SecB_dic, data_dir=data_dir)
print(kf)

#
# fit_dic = yaml.load(Path('batch_fit_config.yaml').read_text())
# fit_dic['epochs'] = 10
# fit_result = do_fitting_from_yaml()