''' convert_data.py 12june2024 LMT
Function to convert data tables exported from HDExaminer to the processed DynamX format pyHDX expects
this will leave any extra columns, but will chop out the MAX time points after processing

Example use case:

# 'all_results.csv' is a summary uptake table exported from HDExaminer
# requires a <protein state>.fasta file for each Protein State 

import pandas as pd
import numpy as np
import os

proj_dir = ''
hdexdf = pd.read_csv(os.path.join(proj_dir,'all_results.csv'))

pepdata = pd.DataFrame()
pepdata = hdexa_to_pyhdx(hdexdf)

hdxm = {}
for mutant in pepdata['state'].unique():
    fasta_sequence =  SeqIO.parse(open(os.path.join(proj_dir,str(mutant)+'.fasta')),'fasta')
    for fasta in fasta_sequence:
        sequence = str(fasta.seq)
        #sequence[mutant] = str(fasta.seq)
    use_data = pepdata.copy()[(pepdata["state"]==mutant) & (pepdata["quality"]!="Low")]
    hdxm[mutant] = HDXMeasurement(use_data, temperature=303.15, pH=8., sequence=sequence)

'''


import pandas as pd
import numpy as np


def hdexa_to_pyhdx(data,d_percentage=0.85,protein='protein'):
    drop_first=2

    def _time_to_sec(tp,tpunit):
        return  tp * np.power(60.0,'smh'.find(tpunit[0])) 
    if '# Deut' in data.columns:
        data = data.rename(columns={"# Deut":"#D"})
        data['#D'] = data['#D'].fillna(0.0)
        data['#D'] = data['#D'].astype(float)
    if 'Deut %' in data.columns:
        data = data.rename(columns={"Deut %":"%D"})
        data['%D'] = data['%D'].fillna(0.0)
        data['%D'] = data['%D'].astype(float)
    if 'Deut Time' in data.columns:
        data.loc[data['Deut Time'] == 'FD','Deut Time'] = '1e6s'
        data['time unit'] = data['Deut Time'].str[-1]
        data['Deut Time (sec)'] = data['Deut Time'].str[:-1].astype(float)
        data['Deut Time (sec)'] = data.apply(lambda x: _time_to_sec(tp=x['Deut Time (sec)'],tpunit=x['time unit']),axis=1)
        data.loc[data['Deut Time (sec)'] == 1e6,'Deut Time (sec)'] = 'MAX'
    if 'Protein' not in data.columns:
        data['Protein'] = protein


    pyhdx_cols = ['start', 'end' ,'stop' ,'sequence', 'state', 'exposure' ,'uptake' ,'maxuptake',
                'fd_uptake' ,'fd_uptake_sd' ,'nd_uptake' ,'nd_uptake_sd' ,'rfu', 'protein',
                'modification', 'fragment', 'mhp' ,'center' ,'center_sd' ,'uptake_sd' ,'rt',
                'rt_sd' ,'rfu_sd' ,'_sequence' ,'_start' ,'_stop' ,'ex_residues',
                'uptake_corrected']
    data = data.rename(columns={
                "Protein State":"state",
                "Protein":"protein",
                "Start":"start",
                "End":"end",
                "Sequence":"_sequence",
                "Peptide Mass":"mhp",
                "RT (min)":"rt",
                "Deut Time (sec)":"exposure",
                "maxD":"maxuptake",
                "Theor Uptake #D":"uptake_corrected",
                "#D":"uptake",
                "%D":"rfu",
                "Conf Interval (#D)":"rfu_sd",
                "#Rep":"rep",
                "Confidence":"quality",
                "Stddev":"center_sd",
                #"p"
    })

    missing = list(set(pyhdx_cols)-set(data.columns))
    for mcol in missing:
        data[mcol] = np.nan
        if mcol == "rfu_sd": data[mcol] = 0.05 #set 5% error as dummy value

    data['rfu']=data['rfu']/100.
    data.loc[data['exposure']=="0",'rfu_sd']=0.0
    data['stop']=data['end']+1
    data['sequence']=data["_sequence"].copy()
    data['sequence']=[s.replace("P", "p") for s in data["sequence"]]
    # Find the total number of n terminal / c_terminal residues to remove from pyhdx/process.py
    n_term = np.array([len(seq) - len(seq[drop_first:].lstrip("p")) for seq in data["sequence"]])
    c_term = np.array([len(seq) - len(seq.rstrip("p")) for seq in data["sequence"]])
    data["sequence"] = ["x" * nt + s[nt:] for nt, s in zip(n_term, data["sequence"])]
    data["_start"] = data["start"] + n_term
    data["_stop"] = data["stop"] - c_term
    ex_residues = (np.array([len(s) - s.count("x") - s.count("p") for s in data["sequence"]])* d_percentage)
    data["ex_residues"] = ex_residues
    data["uptake_sd"]=data["center_sd"]
    data["nd_uptake"]=0.0
    data["nd_uptake_sd"]=0.0
    data["modification"]=float("nan")
    data["fragment"]=float("nan")
    # upeps = data[data["exposure"]=="0"]["_sequence"].unique()
    # fpeps = data[data["exposure"]=="MAX"]["_sequence"].unique()
    # good_peps = np.array(list(set(upeps) & set(fpeps)))
    #peps = data["_sequence"].unique()
    states = data["state"].unique()
    data["fd_uptake"]="novalue"
    data["fd_uptake_sd"]="novalue"

    for state in states:
        peps = data[data["state"]==state]["_sequence"].unique()
        for pep in peps:
            fd_up = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['uptake'].iat[0]
            fd_up_sd = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['center_sd'].iat[0]
            data.loc[data["_sequence"]==pep, "fd_uptake"]=fd_up
            data.loc[data["_sequence"]==pep, "fd_uptake_sd"]=fd_up_sd
    data["center"]=data["mhp"]+data["uptake"]
    data["rt_sd"]=0.05 #dummy value

    data['uptake_corrected_orig'] = data['uptake_corrected'] #sometimes the HDExaminer output value is incorrect
    data['uptake_corrected'] = data["rfu"]*data['maxuptake'] #       so revert to conversion from the rfu

    
    data = data[data["exposure"] != "MAX"]
    data = data[data["fd_uptake"] != 0]
    data = data[~data["uptake"].isna()]
    data["exposure"]=data["exposure"].astype(float)

    new_columns = [col for col in pyhdx_cols if col in data.columns] + [col for col in data.columns if col not in pyhdx_cols]
    return data[new_columns]
