
######################################
# params_loader.py by Sol Samuels
# last edited, 7/2/2026
######################################

# Script to load RBF model parameters from RBF_ModelParameters.xlsx for use by 
# rbf_functions.py . This script must be in the same directory as all other main
# and support scripts/
# This is a support script for RBF_NN.py (head script for RBF-NN implementation)

# If a new state is added to the RBF-NN superposition, the directory "_STATE_ALIAS"
# must be edited to add the new state. The new addition should correspond the python script
# name of the state to the name of the state in RBF_ModelParameters.xlsx
# i.e. "{Python Name for State : Xlsx name for State}"


import pandas as pd
import tensorflow as tf
import keras.backend as K

PARAMS_XLSX = "/users/sol-sam/SolSamuels_AtlasWorkFiles/MANN_Scripts/ExampleRun/RBF_ModelParameters.xlsx"

# Load once when the module is imported
_df = pd.read_excel(PARAMS_XLSX)
_df.columns = [c.strip() for c in _df.columns]

# Map short names to spreadsheet labels
_STATE_ALIAS = {
    "H-1": "H^{-1}",
    "H+1": "H^{+1}",
    "Li0": "Li^{0}",
    "Li+1": "Li^{+1}",
    "Li-1": "Li^{-1}",
    "Liexc1": "Li^{exc1}",
    "Liexc2": "Li^{exc2}",
    "Liexc3": "Li^{exc3}",
    "F0": "F^{0}",
    "F+1": "F^{+1}",
    "F-1": "F^{-1}",
    "F-2": "F^{-2}",
    "Fexc1": "F^{exc1}",
    "Fexc2": "F^{exc2}",
    "C0": "C^{0}",
    "C+1": "C^{+1}",
    "C-1": "C^{-1}",
    "C+2": "C^{+2}",
    "Cexc1": "C^{exc1}",
    "Cexc2": "C^{exc2}",
    "Cexc3": "C^{exc3}",
    "O0": "O^{0}",
    "O+1": "O^{+1}",
    "O-1": "O^{-1}",
    "O-2": "O^{-2}",
    "Oexc1": "O^{exc1}",
    "Oexc2": "O^{exc2}",
    "N0": "N^{0}",
    "N+1": "N^{+1}",
    "N-1": "N^{-1}",
}

_REQUIRED_KEYS = [
    "Z","A0","B0","C0","D0","E0",
    "alpha","beta","bb","mm","gamma","eta","delta","pp","qq","tt"
]

def params(state: str, dtype=None):
    if dtype is None:
        dtype = tf.as_dtype(K.floatx())  # typically float32
    key = _STATE_ALIAS.get(state, state)
    rows = _df.loc[_df["Atomic state"] == key]
    if rows.empty:
        raise KeyError(f"State '{state}' not found (looked for '{key}').")
    row = rows.iloc[0].fillna(0.0)
    param_dict = {k: tf.constant(float(row[k]), dtype=dtype) for k in _REQUIRED_KEYS}
    return param_dict
