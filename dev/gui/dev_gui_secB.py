"""
Reload SecB and fitted data and launch  GUI
Run local_cluster.py in anothor thread

"""

from pyhdx.fileIO import read_dynamx, csv_to_dataframe
from pyhdx import PeptideMasterTable
import pickle
from pyhdx.web.apps import main_app
from pyhdx.web.base import DEFAULT_COLORS, STATIC_DIR
from pyhdx.web.utils import load_state
from pyhdx.web.sources import DataSource
from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_protein
import panel as pn
import numpy as np
from pathlib import Path
import pandas as pd
import yaml

ctrl = main_app()
directory = Path(__file__).parent
root_dir = directory.parent.parent
data_dir = root_dir / 'tests' / 'test_data' / 'input'
test_dir = directory / 'test_data'

fpath_1 = root_dir / 'tests' / 'test_data' / 'ecSecB_apo.csv'
fpath_2 = root_dir / 'tests' / 'test_data' / 'ecSecB_dimer.csv'

yaml_dict = yaml.safe_load(Path(data_dir / 'data_states.yaml').read_text())

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


def reload_dashboard():
    data_objs = {k: yaml_to_hdxm(v, data_dir=data_dir) for k, v in yaml_dict.items()}
    for k, v in data_objs.items():
        v.metadata['name'] = k
    ctrl.data_objects = data_objs

    source = ctrl.sources['dataframe']
    for ds in ['peptides', 'peptides_mse', 'd_calc', 'rfu', 'rates', 'global_fit', 'losses']:
        df = csv_to_protein(test_dir / f'{ds}.csv')
        source.add_df(df, ds)

    #Temporary workaround for comment characters in csv files
    ds = 'colors'
    df = pd.read_csv(test_dir / f'{ds}.csv', header=[0, 1, 2], index_col=0,
                     skiprows=3)
    source.add_df(df, ds)

    # rates = csv_to_protein(test_dir / 'rates.csv').df
    #
    # fit = csv_to_protein(test_dir / 'global_fit.csv').df
    # colors = csv_to_protein(test_dir / 'colors.txt').df
    # peptides = csv_to_dataframe(test_dir / 'peptides.txt')
    #
    # source = ctrl.sources['dataframe']
    # source.add_df(rates, 'rates')
    # source.add_df(peptides, 'peptides')
    # source.add_df(fit, 'global_fit')
    # #source.add_df(colors, 'colors')

    #ctrl.sources['dataframe'].updated = True

    # fit_control = ctrl.control_panels['FitControl']
    # fit_control.epochs = 100
    # fit_control.fit_mode = 'Single'
    # fit_control.fit_name = 'new_global_fit_test_123'
    #
    # ngl = ctrl.views['protein']
    # ngl.ngl_view.pdb_string = Path(test_dir / '1qyn.pdb').read_text()


def init_dashboard():
    for k, v in yaml_dict.items():
        load_state(ctrl, v, data_dir=data_dir, name=k)

    # file_input = ctrl.control_panels['PeptideFileInputControl']
    # file_input.input_files = files
    # file_input.fd_state = 'Full deuteration control'
    # file_input.fd_exposure = 0.167*60
    # file_input.pH = 8
    # file_input.temperature = 273.15 + 30
    # file_input.d_percentage = 90.
    #
    # file_input.exp_state = 'SecB WT apo'
    # file_input.dataset_name = 'SecB_tetramer'
    # file_input._action_add_dataset()
    #
    # file_input.exp_state = 'SecB his dimer apo'
    # file_input.dataset_name = 'SecB_dimer'  # todo catch error duplicate name
    # file_input._action_add_dataset()

    # initial_guess = ctrl.control_panels['InitialGuessControl']
    # initial_guess._action_fit()
#
    fit_control = ctrl.control_panels['FitControl']
    fit_control.epochs = 10

    fit_control.r1 = 0.05
    fit_control.r2 = 0.1
    fit_control.epochs = 200000
    fit_control.stop_loss = 0.001
    fit_control.patience = 100
    fit_control.learning_rate = 100


    ngl = ctrl.views['protein']
    ngl.ngl_view.pdb_string = Path(test_dir / '1qyn.pdb').read_text()

#     fit_control.fit_mode = 'Batch'
#     fit_control._action_fit()
#     fit_control._do_fitting()
# #
#     classification = ctrl.control_panels['ClassificationControl']
#     classification.widgets['select_1'].value = '*'
#     classification.widgets['select_2'].value = 'deltaG'
#
#     classification.mode = 'Continuous'
#     classification._action_linear()
#     classification.color_set_name = 'colorset test'
#     classification._action_add_colorset()
# #
#
# file_export = ctrl.control_panels['FileExportControl']
# sio = file_export.table_export_callback()


#if __name__ == '__main__':
#pn.state.onload(reload_dashboard)
pn.state.onload(init_dashboard)

pn.serve(ctrl.template, show=True
         , static_dirs={'pyhdx': STATIC_DIR})

#ctrl.template.servable()
