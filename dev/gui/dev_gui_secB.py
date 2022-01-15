"""
Reload SecB and fitted data and launch  GUI
Run local_cluster.py in anothor thread

"""

import sys
from pathlib import Path

import pandas as pd
import panel as pn
import yaml

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.fileIO import csv_to_protein
from pyhdx.web.apps import main_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.utils import load_state, fix_multiindex_dtypes

sys._excepthook = sys.excepthook

import traceback as tb
def my_exception_hook(exctype, value, traceback):
    # Print the error and traceback
    # https://stackoverflow.com/questions/43039048/pyqt5-fails-with-cryptic-message/43039363#43039363
    tb.print_tb(traceback, file=sys.stdout)
    print(exctype, value, traceback)

    tb.print_stack()
    print(traceback.format_exc())
    # or
    print(sys.exc_info()[2])
    # Call the normal Exception hook after
    sys._excepthook(exctype, value, traceback)
    sys.exit(1)

# Set the exception hook to our wrapping function
sys.excepthook = my_exception_hook


ctrl, tmpl = main_app()


directory = Path(__file__).parent.resolve()
root_dir = directory.parent.parent
data_dir = root_dir / 'tests' / 'test_data' / 'input'
test_dir = directory / 'test_data'

fpath_1 = root_dir / 'tests' / 'test_data' / 'ecSecB_apo.csv'
fpath_2 = root_dir / 'tests' / 'test_data' / 'ecSecB_dimer.csv'
fitresult_dir = root_dir / 'tests' / 'test_data' / 'output' / 'ecsecb_tetramer_dimer'

yaml_dict = yaml.safe_load(Path(data_dir / 'data_states.yaml').read_text())
pdb_string = (test_dir / '1qyn.pdb').read_text()


def reload_tables():
    src = ctrl.sources['main']
    src.tables['peptides'] = csv_to_dataframe(test_dir / 'peptides.csv')
    src.tables['rfu_residues'] = csv_to_dataframe(test_dir / 'rfu_residues.csv')
    src.tables['dG_fits'] = csv_to_dataframe(test_dir / 'dG_fits.csv')
    src.tables['ddG_comparison'] = csv_to_dataframe(test_dir / 'ddG_comparison.csv')
    src.tables['rates'] = csv_to_dataframe(test_dir / 'rates.csv')
    src.param.trigger('updated')


    ctrl.views['protein'].object = pdb_string

def reload_dashboard():
    data_objs = {k: yaml_to_hdxm(v, data_dir=data_dir) for k, v in yaml_dict.items()}
    for k, v in data_objs.items():
        v.metadata['name'] = k

    source = ctrl.sources['dataframe']
    for ds in ['peptides', 'peptides_mse', 'd_calc', 'rfu', 'rates', 'global_fit', 'losses']:
        df = csv_to_protein(test_dir / f'{ds}.csv')
        source.add_df(df, ds)

    #Temporary workaround for comment characters in csv files
    ds = 'colors'
    df = pd.read_csv(test_dir / f'{ds}.csv', header=[0, 1, 2], index_col=0,
                     skiprows=3)
    source.add_df(df, ds)

def init_dashboard():
    for k, v in yaml_dict.items():
        load_state(ctrl, v, data_dir=data_dir, name=k)

    src = ctrl.sources['main']
    fit_control = ctrl.control_panels['FitControl']
    fit_control.epochs = 10

    fit_control.r1 = 0.05
    fit_control.r2 = 0.1
    fit_control.epochs = 200000
    fit_control.stop_loss = 0.001
    fit_control.patience = 100
    fit_control.learning_rate = 100

    # ngl = ctrl.views['protein']
    # ngl._ngl.pdb_string = Path(test_dir / '1qyn.pdb').read_text()
    # ctrl.views['protein'].object = pdb_string
    #
    # pdbe = ctrl.views['protein']

    fit_result = load_fitresult(fitresult_dir)
    src.add(fit_result, 'fit_1')

    diff = ctrl.control_panels['DifferentialControl']
    diff._action_add_comparison()


#if __name__ == '__main__':
#pn.state.onload(reload_dashboard)
#pn.state.onload(reload_tables)
pn.state.onload(init_dashboard)


if __name__ == '__main__':
    pn.serve(tmpl, show=True, static_dirs={'pyhdx': STATIC_DIR})

elif __name__.startswith('bokeh_app'):
    tmpl.servable()
