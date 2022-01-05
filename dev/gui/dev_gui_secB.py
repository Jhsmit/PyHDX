"""
Reload SecB and fitted data and launch  GUI
Run local_cluster.py in anothor thread

"""

import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import pandas as pd
import panel as pn
import yaml

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.fileIO import csv_to_protein
from pyhdx.local_cluster import default_client
from pyhdx.web.apps import main_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.cache import MemoryCache
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

cwd = Path(__file__).parent.resolve()
root_dir = cwd.parent.parent
test_dir = cwd / 'test_data'

fpath_1 = root_dir / 'tests' / 'test_data' / 'ecSecB_apo.csv'
fpath_2 = root_dir / 'tests' / 'test_data' / 'ecSecB_dimer.csv'
fitresult_dir = root_dir / 'tests' / 'test_data' / 'output' / 'ecsecb_tetramer_dimer'


data_dir = root_dir / 'tests' / 'test_data' / 'input'
#data_dir = cwd / 'rinky_data'

yaml_dict = yaml.safe_load(Path(data_dir / 'data_states.yaml').read_text())
pdb_string = (test_dir / '1qyn.pdb').read_text()

# fpaths = [fpath_1, fpath_2]
# files = [p.read_bytes() for p in fpaths]
#
#
# d1 = {
#     'filenames': ['ecSecB_apo.csv', 'ecSecB_dimer.csv'],
#     'd_percentage': 95,
#     'control': ('Full deuteration control', 0.167),
#     'series_name': 'SecB WT apo',
#     'temperature': 30,
#     'temperature_unit': 'celsius',
#     'pH': 8.,
#     'c_term': 165
# }
#
# d2 = {
#     'filenames': ['ecSecB_apo.csv', 'ecSecB_dimer.csv'],
#     'd_percentage': 95,
#     'control': ('Full deuteration control', 0.167),
#     'series_name': 'SecB his dimer apo',
#     'temperature': 30,
#     'temperature_unit': 'celsius',
#     'pH': 8.,
#     'c_term': 165
# }

#yaml_dicts = {'testname_123': d1, 'SecB his dimer apo': d2}


def reload_tables():
    test_dir = cwd / 'test_data'
    src = ctrl.sources['main']

    df = csv_to_dataframe(test_dir / 'peptides.csv')
    # names = df.columns.names
    # df = df.convert_dtypes()
    # df.columns.names = names
    table_names = [
            'rfu_residues.csv',
            'rates.csv',
            'peptides.csv',
            'dG_fits.csv',
            'ddG_comparison.csv',
            'd_calc.csv',
            'loss.csv',
            'peptide_mse.csv'
    ]
    for name in table_names:
        try:
            df = csv_to_dataframe(test_dir / name)
            df.columns = fix_multiindex_dtypes(df.columns)
            src.tables[name.split('.')[0]] = df
        except Exception as e:
            print(e)
            print('not loaded:', name)


    # src.tables['peptides'] = fix_index_dtypes(df)  #csv_to_dataframe(test_dir / 'peptides.csv')
    #
    # src.tables['rfu_residues'] = csv_to_dataframe(test_dir / 'rfu_residues.csv')
    # src.tables['dG_fits'] = csv_to_dataframe(test_dir / 'dG_fits.csv')
    # src.tables['ddG_comparison'] = csv_to_dataframe(test_dir / 'ddG_comparison.csv')
    # src.tables['peptide_mse'] = csv_to_dataframe(test_dir / 'peptide_mse.csv')
    # src.tables['d_calc'] = csv_to_dataframe(test_dir / 'd_calc.csv')
    # #src.tables['rates'] = csv_to_dataframe(test_dir / 'rates.csv')
    src.param.trigger('updated')


    #ctrl.views['protein'].object = pdb_string

#ctrl.views['protein'].object = pdb_string


def reload_dashboard():
    data_objs = {k: yaml_to_hdxm(v, data_dir=data_dir) for k, v in yaml_dict.items()}
    for k, v in data_objs.items():
        v.metadata['name'] = k

    source = ctrl.sources['main']
    for ds in ['peptides', 'peptides_mse', 'd_calc', 'rfu', 'rates', 'global_fit', 'losses']:
        df = csv_to_protein(test_dir / f'{ds}.csv')
        source.add_df(df, ds)

    #Temporary workaround for comment characters in csv files
    ds = 'colors'
    df = pd.read_csv(test_dir / f'{ds}.csv', header=[0, 1, 2], index_col=0,
                     skiprows=3)
    source.add_df(df, ds)


def init_dashboard(ctrl):
    for k, v in yaml_dict.items():
        load_state(ctrl, v, data_dir=data_dir, name=k)

    # k = next(iter(yaml_dict.keys()))
    # load_state(ctrl, yaml_dict[k], data_dir=data_dir, name=k)

    src = ctrl.sources['main']
    fit_control = ctrl.control_panels['FitControl']

    fit_control.r1 = 0.05
    fit_control.r2 = 0.1
    fit_control.epochs = 200
    fit_control.patience = 100

    # ngl = ctrl.views['protein']
    # ngl._ngl.pdb_string = Path(test_dir / '1qyn.pdb').read_text()
    fit_result = load_fitresult(fitresult_dir)

    src.add(fit_result, 'gibbs_fit_1')

    pdb_src = ctrl.sources['pdb']
    pdb_src.add_from_string(pdb_string, '1qyn')

    diff = ctrl.control_panels['DifferentialControl']
    #diff._action_add_comparison()

    # ctrl.views['protein'].object = pdb_string
    #
    # fit_result = load_fitresult(fitresult_dir)
    # src.add(fit_result, 'fit_1')


#pn.state.onload(reload_dashboard)
#pn.state.onload(reload_tables)

if __name__ == '__main__' or __name__.startswith('bokeh_app'):
    # sys._excepthook = sys.excepthook
    #
    # import traceback as tb
    # def my_exception_hook(exctype, value, traceback):
    #     # Print the error and traceback
    #     # https://stackoverflow.com/questions/43039048/pyqt5-fails-with-cryptic-message/43039363#43039363
    #     tb.print_tb(traceback, file=sys.stdout)
    #     print(exctype, value, traceback)
    #
    #     tb.print_stack()
    #     print(traceback.format_exc())
    #     # or
    #     print(sys.exc_info()[2])
    #     # Call the normal Exception hook after
    #     sys._excepthook(exctype, value, traceback)
    #     sys.exit(1)
    #
    #
    #
    # # Set the exception hook to our wrapping function
    # sys.excepthook = my_exception_hook

    #init_dashboard()
    print('starting')
    #executor = default_client(asynchronous=True)
    executor = ProcessPoolExecutor(max_workers=2) # doesnt work / extremely slow (however one-time, torch works fine)
    cache = MemoryCache(max_items=2000)

    ctrl, tmpl = main_app(executor, cache)

    loader = init_dashboard
    func = partial(loader, ctrl)

    pn.state.onload(func)
    if __name__ == '__main__':
        pn.serve(tmpl, show=True, static_dirs={'pyhdx': STATIC_DIR})
    elif  __name__.startswith('bokeh_app'):
        tmpl.servable()
        # panel serve --show --autoreload --static-dirs pyhdx=C:\Users\jhsmi\pp\PyHDX\pyhdx\web\static

#ctrl.template.servable()


