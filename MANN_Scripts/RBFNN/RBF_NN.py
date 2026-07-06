"""
Versatile Radial Basis Function Neural Network (RBF-NN) for Diatomics

Made by Sol Samuels
Last edited: 7/6/2026

This script implements a TensorFlow, custom Keras layer, radial basis function neural network 
designed to determine weights for an ensemble atom in molecule representation [1] 
for molecular electron densities for diatomic systems as a function of internuclear
separation (R). The model is trained on precomputed molecular electron density grids (Gaussian or
similar electronic structure outputs) and determines a weighted atomic-state decomposition that
reconstructs the full molecular density.


[1]  Susan R Atlas. Embedding quantum statistical excitations in a classical force field.
The Journal of Physical Chemistry A, 125(17):3760-3775, 2021
"""

# ---------------------------------
# Imports Required:
# ---------------------------------
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad
import keras.backend as K
import sys
from scipy.integrate import trapz
from scipy.interpolate import RegularGridInterpolator
import os, datetime
import re
from datetime import date
from datetime import datetime
import time

from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import Callback


# The following 3 imports require rbf_functions.py and params_loader.py to be in the job directory
# The file RBF_ModelParameters.xlsx must also be in the job directory
from rbf_functions import M5, M6, M7
from rbf_functions import hydrogen_density as HDen
from params_loader import params


print('imports complete')

total_start = time.perf_counter()

TRIAL_ENV = os.environ.get("TRIAL", "default")
trial = TRIAL_ENV



# ---------------------------------
# SETUP OF DIATOMIC SYSTEM:
# ---------------------------------
#    Change to match the diatomic system being studied and the level of theory.
#    Note that these settings must match pre-calculated training data calculated 
#    using Gaussian (molecular electron density for various internuclear separations R)


level_of_theory = "MP4"     # Level of theory of training data
basis_set = "aug-cc-pVQZ"   # Basis set of training data
molecule = "CO"             # Molecule being studied
AtomA = "C"                 # Atom at origin
AtomB = "O"                 # Atom at z=R


# ---------------------------------
# SETUP OF RBF_NN PARAMETERS:
# ---------------------------------
#    These parameters may be changed to affect learning of the NN


# Training Parameters
lr = 0.01                   # learning-rate input into ADAM
coarseness = "Coarse"       # Training grid level of coarseness (Gaussian determined)
batch_s = 1024              # Training batch size (Try Coarse=1024, Medium=8192, Fine=65536)
start_dcy = 0.0001          # L2 weight decay parameter (don't go too large)
end_dcy = start_dcy         # Change if want different L2 weight decay parameter in short-range (small R)
epoch_num = 30000           # Maximum number of epochs before training is ended
stopping_patience = 500     # Number of epochs with no improvement in loss before training is ended
epsilon = 1e-2              # Tolerance for weight decay switch (ignore if end_dcy = start_dcy)

# Loss function parameters:
loss_fitparam = 2e4         # weight for mean(square(rho_true - rho_pred)) term
loss_chargeparam = 100      # weight for total charge constraint
loss_AtomAsum = 100         # weight for sum of atom A weights = 1 constraint
loss_AtomBsum = 100         # weight for sum of atom B weights = 1 constraint


# Fraction of the values of R included in the weight calculation (0 to 1)
    # if set to 1, all available values of R (those with molecular density grids calculated) will be
    # included in the RBF-NN weight determination. If set to 0.5, every other available distance will be skipped, etc.
fraction_sep_dict = 1.0  # Reduce if want a quick 'test' run


# Define requested states: (name, charge, model function)
ATOM_A_STATES = [
    {"name": "0",    "charge":  0, "func": M5},
    {"name": "+1",   "charge":  1, "func": M5},
    {"name": "-1",   "charge": -1, "func": M6},
    {"name": "+2",   "charge":  2, "func": M5},
    {"name": "exc1", "charge":  0, "func": M7},
    {"name": "exc2", "charge":  0, "func": M7},
    {"name": "exc3", "charge":  0, "func": M7},
]

ATOM_B_STATES = [
    {"name": "0",    "charge":  0, "func": M5},
    {"name": "+1",   "charge":  1, "func": M5},
    {"name": "-1",   "charge": -1, "func": M6},
    {"name": "-2",   "charge": -2, "func": M6},
    {"name": "exc1", "charge":  0, "func": M7},
    {"name": "exc2", "charge":  0, "func": M7},
]

# Example for Hydrogen states:
#ATOM_A_STATES = [
 #   {"name": "1s", "charge": 0, "func": HDen,
 #    "kwargs": {"n": 1, "l": 0}},
 #   {"name": "2s", "charge": 0, "func": HDen,
 #    "kwargs": {"n": 2, "l": 0}},
 #   {"name": "2p", "charge": 0, "func": HDen,
 #    "kwargs": {"n": 2, "l": 1}},
#]


N_A = len(ATOM_A_STATES)
N_B = len(ATOM_B_STATES)
NODES = N_A + N_B           # Number of nodes = number of states


# ---------------------------------
# CREATION OF OUTPUT DIRECTORY
# ---------------------------------

# Each run of the RBF-NN results in outputs in a new directory
TASK_NAME = os.environ.get("TASK_NAME", f"{molecule}_RBF")
RUN_TAG = os.environ.get("RUN_TAG", None)

OUTDIR = os.environ.get(
    "OUTDIR",
    os.path.join(
        os.getcwd(),
        RUN_TAG if RUN_TAG else f"{TASK_NAME}_T{TRIAL_ENV}_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')}"
    ),
)
os.makedirs(OUTDIR, exist_ok=True)
os.chdir(OUTDIR)

# ---------------------------------
# PRINT OUT RBF_NN PARAMETERS:
# ---------------------------------

print("\nRBF-NN Initializing:")
print("\tMolecular System: ", molecule)
print(f"\tLevel of Theory / Basis Set: {level_of_theory}/{basis_set} ")
print("\tTrial: ", trial)
print("\tTraining Grid Coarseness: ", coarseness)
print(f"\tWriting outputs to: {OUTDIR}")


print("\nParameter Information:")
print("\t Learning rate: ", lr)
print("\t Batch Size: ", batch_s)
print("\t Long-range Weight Decay: ", start_dcy)
print("\t Short-range Weight Decay: ", end_dcy)
print("\t Weight Dcy switch eps. tolerance (ignore if end_dcy = start_dcy): ", epsilon)
print("\t Max Epoch number per value of R: ", epoch_num)
print("\t Epoch Stopping Patience: ", stopping_patience)
print("\t Loss Function RSME weight: ", loss_fitparam)
print("\t Loss Function charge constraint weight: ", loss_chargeparam)
print("\t Loss Function Atom A sum constraint: ", loss_AtomAsum)
print("\t Loss Function Atom B sum constraint: ", loss_AtomBsum)





#File path for density files:
if coarseness == "Coarse":
    main_path = f"/carc/scratch/projects/susie/susie2016319-3tb/GaussianMolecularData/{molecule}/{level_of_theory}/{coarseness}/{molecule}_{level_of_theory}_{basis_set}_seq/"
if coarseness == "Fine":
    main_path = f"/easley/scratch/users/sol-sam/FineMolecularDensities/{molecule}/{level_of_theory}/{coarseness}/{molecule}_{level_of_theory}_{basis_set}_seq/"




# Internuclear separations (in Angstrom)
# Build separation list automatically from directory naming:
sep_dict = []

for entry in os.listdir(main_path):
    if os.path.isdir(os.path.join(main_path, entry)) and entry.endswith("A"):
        try:
            R = float(entry[:-1])   # strip trailing "A"
            sep_dict.append((R, entry[:-1]))
        except ValueError:
            pass

# sort largest -> smallest (large R to small R)
sep_dict.sort(reverse=True)

# keep only the string values used elsewhere in the script
sep_dict = [s for _, s in sep_dict]
            
# --------------------------
# Include only a percentage of R values if requested
# --------------------------

# Keep an evenly distributed fraction
if fraction_sep_dict < 1.0:
    n = len(sep_dict)
    keep = max(1, round(n * fraction_sep_dict))

    indices = np.unique(np.round(np.linspace(0, n - 1, keep)).astype(int))
    sep_dict = [sep_dict[i] for i in indices]
    print(f"\n{fraction_sep_dict}% of values of R will be included")
else:
    print("\nFraction of internuclear separations requested is not less than 1. All possible values of R are included.")

print(f"Number of values of R to be evaluated: {len(sep_dict)}")

# Initial guess for weights
# Default initial guess is both atoms are neutral
ini_weights = np.zeros(NODES)
ini_weights[0] = 1
ini_weights[N_A] = 1

# Header for mas_data excel file output which excludes results
mas_header = [
    'R (A)',
    'R (a.u.)',
    'batch size',
    'lr',
    'dcy const'
]

mas_header += [
    f"{AtomA}_{state['name']}"
    for state in ATOM_A_STATES
]

mas_header += [
    f"{AtomB}_{state['name']}"
    for state in ATOM_B_STATES
]

mas_header += [
    f'{AtomA} Eff.',
    f'{AtomB} Eff.',
    'Tot. Ch',
    'RMSE',
    'epoch len.'
]

mas_data = [mas_header]



# ---------------------------------
# DEFINING THE RBF_NN IN KERAS:
# ---------------------------------
#     This Class format defines the custom RBF-NN layer
#     Should be robust between changes of the diatomic system being studied, as long as
#     all requested states have been fit with parameters listed in RBF_ModelParameters.xlsx

class RadialBasisFunction(tf.keras.layers.Layer):
    def __init__(self, nodes, R, weights):
        super(RadialBasisFunction, self).__init__()
        self.units = nodes
        self.R = R
        
        self.epsilon = 1e-8  # A small number
        
        self.initial_weights = weights

        self.A_states = []
        self.B_states = []

        for state in ATOM_A_STATES:
            entry = {
                "name": state["name"],
                "charge": state["charge"],
                "func": state["func"]
            }

            if state["func"] == HDen:
                entry["kwargs"] = state["kwargs"]
            else:
                entry["params"] = params(f"{AtomA}{state['name']}")
        
            self.A_states.append(entry)

        for state in ATOM_B_STATES:
            entry = {
                "name": state["name"],
                "charge": state["charge"],
                "func": state["func"]
            }

            if state["func"] == HDen:
                entry["kwargs"] = state["kwargs"]
            else:
                entry["params"] = params(f"{AtomB}{state['name']}")
        
            self.B_states.append(entry)

        
    def build(self, input_shape):
        # Adds a weight variable to the layer
        self.A_weights = []
        for i in range(N_A):
            self.A_weights.append(
                self.add_weight(
                    name=f"wA{i}",
                    shape=(1,),
                    initializer=tf.keras.initializers.Constant(
                        value=self.initial_weights[i]
                    ),
                    trainable=True,
                    constraint=tf.keras.constraints.NonNeg()
                )
            )

        self.B_weights = []
        for i in range(N_B):
            self.B_weights.append(
                self.add_weight(
                    name=f"wB{i}",
                    shape=(1,),
                    initializer=tf.keras.initializers.Constant(
                        value=self.initial_weights[N_A+i]
                    ),
                    trainable=True,
                    constraint=tf.keras.constraints.NonNeg()
                )
            )

        super().build(input_shape) 
        
        
    def call(self, inputs):
        coor = inputs
        RBF_0 = 0.0
        RBF_1 = 0.0

        for w, state in zip(self.A_weights, self.A_states):
            if state["func"] == HDen:
                RBF_0 += w * HDen(
                    self,
                    coor,
                    0,
                    state["kwargs"]["n"],
                    state["kwargs"]["l"]
                )
            else:
                RBF_0 += w * state["func"](
                    self,
                    coor,
                    0,
                    state["params"]
                )

        for w, state in zip(self.B_weights, self.B_states):
            if state["func"] == HDen:
                RBF_1 += w * HDen(
                    self,
                    coor,
                    self.R,
                    state["kwargs"]["n"],
                    state["kwargs"]["l"]
                )
            else:
                RBF_1 += w * state["func"](
                    self,
                    coor,
                    self.R,
                    state["params"]
                )
        
        RBF_r = RBF_0 + RBF_1
        
        return RBF_r
    



# ---------------------------------
# EARLY STOPPING ROUTINE
# ---------------------------------
#     Monitors loss and ends training if there is no improvement after a certain patience
#     Returns the best weights with the lowest loss

class EarlyStoppingWithBestWeights(Callback):
    def __init__(self, monitor='val_loss', min_delta=0, patience=5, n_average=20, factor=0.9, verbose=0, mode='auto', baseline=None, start_epoch=0):
        super(EarlyStoppingWithBestWeights, self).__init__()

        self.monitor = monitor
        self.min_delta = min_delta
        self.patience = patience
        self.n_average = n_average
        self.factor = factor
        self.verbose = verbose
        self.wait = 0
        self.stopped_epoch = 0
        self.best = None
        self.best_weights = None  # New attribute to store the best weights
        self.mode = mode
        self.baseline = baseline
        self.loss_history = []
        self.start_epoch = start_epoch

    def on_train_begin(self, logs=None):
        self.best = float('inf')

    def on_epoch_end(self, epoch, logs=None):
        if epoch < self.start_epoch:
            return  # Do not evaluate until start_epoch is reached

        current = logs.get(self.monitor)
        
        if current is not None:
            if self.mode == 'auto':
                if 'acc' in self.monitor or 'f1' in self.monitor:
                    mode = 'max'
                else:
                    mode = 'min'
            else:
                mode = self.mode
        # Update the loss history and check for improvement
        if current is not None:
            self.loss_history.append(current)
            if current < 5:
                # average_loss = np.mean(self.loss_history[-self.n_average:])
                if current < self.best:
                    self.best = current
                    self.best_weights = self.model.get_weights()  # Update best weights
                    self.wait = 0
                else:
                    self.wait += 1

                if self.wait >= self.patience:
                    self.stopped_epoch = epoch
                    self.model.set_weights(self.best_weights)  # Set model weights to the best weights
                    self.model.stop_training = True
                    if self.verbose > 0:
                        print(f'\nEpoch {epoch + 1}: Early stopping due to no improvement in {self.monitor}. Best weights restored with RMSE {self.best}.')



# ---------------------------------
# CALL WEIGHTS AND CHARGES FUNCTIONS:
# ---------------------------------
#       These functions pull the weights from the keras layer
#       for atom A and atom B, as well as determine the effective charge
#       and sum of the weights for each value of R

def get_layer_weights():

    try:
        weights = model.layers[0].weights
    except NameError:
        weights = model_val.layers[0].weights

    return weights


def split_weights(weights):

    A_weights = weights[:N_A]
    B_weights = weights[N_A:]

    return A_weights, B_weights


def effective_charge(weight_list, state_list):

    return tf.add_n([
        w * state["charge"]
        for w, state in zip(weight_list, state_list)
    ])


def weight_sum(weight_list):

    return tf.add_n(weight_list)




# ---------------------------------
# LOSS FUNCTION:
# ---------------------------------
#     The following function designates the loss function which is called on by the RBF-NN 
#     during training

def loss_fit(rho_true, rho_pred):

    weights = get_layer_weights()

    A_weights, B_weights = split_weights(weights)

    A_charge = effective_charge(A_weights, ATOM_A_STATES)
    B_charge = effective_charge(B_weights, ATOM_B_STATES)

    reduced_loss = (
        loss_fitparam * tf.reduce_mean(tf.square(rho_true - rho_pred))   # MSE evaluation
        + loss_chargeparam * tf.abs(A_charge + B_charge)                 # Total charge constraint
        + loss_AtomAsum * tf.abs(weight_sum(A_weights) - 1)              # weights summing to one constraint
        + loss_AtomBsum * tf.abs(weight_sum(B_weights) - 1)              # for both atom A and B
    )

    return reduced_loss




# ---------------------------------
# TRAINING MONITORING FUNCTIONS:
# ---------------------------------
#     The following functions are used to monitor training success


# Function to calulate RMSE between true and predicted density
def rmse(rho_true, rho_pred):
    return (tf.reduce_mean(tf.square(rho_true - rho_pred)))


# Function to calulate the total charge
def total_charge(rho_true, rho_pred):

    weights = get_layer_weights()

    A_weights, B_weights = split_weights(weights)

    return (
        effective_charge(A_weights, ATOM_A_STATES)
        + effective_charge(B_weights, ATOM_B_STATES)
    )


# Function to calulate the sum of weights for Atom A
def AtomA_sum(rho_true, rho_pred):

    weights = get_layer_weights()

    A_weights, B_weights = split_weights(weights)

    return weight_sum(A_weights)



# Function to calulate the sum of weights for Atom B
def AtomB_sum(rho_true, rho_pred):

    weights = get_layer_weights()

    A_weights, B_weights = split_weights(weights)

    return weight_sum(B_weights)



# Function to calulate the charge on Atom A
def q(rho_true, rho_pred):

    weights = get_layer_weights()

    A_weights, B_weights = split_weights(weights)

    return effective_charge(
        A_weights,
        ATOM_A_STATES
    )



# Function to get value of shift along the z-axis for placement of atom A at the origin
def get_shift(path):
    with open(path, "r") as f:
        text = f.read()

    s = f"with {AtomA} at z="
    i = text.find(s)
    if i == -1:
        return None
    i += len(s)
    j = i
    while j < len(text) and (text[j].isdigit() or text[j] in "+-.eE"):
        j += 1
    return float(text[i:j])



# ---------------------------------
# CALLING THE RBF-NN:
# ---------------------------------
#     

wA_hist = [[] for x in range(N_A)]
wB_hist = [[] for x in range(N_B)]

dipole_pred_lst = []
dipole_train_lst = []


best_loss = float('inf')
dcy = start_dcy
charge_track = []
trans_flag = False


# Iterates though each value of R (from large R to small R)
for key_str in sep_dict:
    sep_start = time.perf_counter()
    key = float(key_str)
    print("\n-------------------------------------------------------")
    print("Separation: ", key_str)
    print("\n")
    #print(charge_track)
    
    
    # Monitors running average for change in weight decay
    if len(charge_track) > 0:
        run_ave = np.mean(charge_track)
        print("Running Ave: ", run_ave)
        
    if len(charge_track) >= 5 and trans_flag == False:
        window = charge_track[-3:]
        ave_ch = np.mean(window)
        print("Average Charge: ", ave_ch)
        print("Reference Charge: ", run_ave)
        if ave_ch - epsilon > run_ave:
            dcy = end_dcy
            print("Decay Changed to :", dcy)
            trans_flag = True
            
        else:
            print("Decay Unchanged")
            
        print('\n')
            
    try:
        epoch_data = [['epoch', 'training_loss', 'q_per_epoch', 'rmse_ep']]

        R = key * 1.8897259886
        
        # README needs to be accessable at this path to get shift along z-axis
        shift_path = f"{key_str}A/README.txt"
        path = main_path + shift_path
        shift = get_shift(path)
        print("Shift: ", shift)

        data = []
        val_data = []
        
        
        # Collects training data from xyz .out density file and applies shift
        with open(main_path + f"{key_str}A/{molecule}_{key_str}A_Den{coarseness}.out", "r") as file:
            for line in file:
                # Split the line by the delimiter (e.g., tab "\t")
                values = line.strip().split(" ")
                realvals = []
                for val in values:
                    if val != '':
                        realvals.append(val)

                # Convert the relevant columns to float
                column1 = float(realvals[0])
                column2 = float(realvals[1]) 
                column3 = float(realvals[2]) - shift   # Applies shift to place atom A at origin
                column4 = float(realvals[3]) 

                #if column4 > 1e-7:
                    #val_data.append([column1, column2, column3, column4])

                data.append([column1, column2, column3, column4])




        Molecule_cube = pd.DataFrame(data, columns=["col1", "col2", "col3","col4"])
        Molecule_cube_arr = Molecule_cube.to_numpy()
        #np.random.shuffle(Molecule_cube_arr)

        Molecule_cube_coor = Molecule_cube_arr[:,:3]
        Molecule_cube_rho = Molecule_cube_arr[:,3]

        coor_train = Molecule_cube_coor
        rho_train = Molecule_cube_rho

        coor_test = Molecule_cube_coor
        rho_test = Molecule_cube_rho

        
        

        # Number of Nodes = Number of atomic states
        nodes = N_A + N_B

        # Defining ADAM
        adam = tf.keras.optimizers.legacy.Adam(learning_rate=lr, decay=dcy, amsgrad=True) 
         
        # Creating a sequential model with RBF layer
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.InputLayer(input_shape=(3,)))
        model.add(RadialBasisFunction(nodes, R, ini_weights)) #will need to add other parameters

        # Compiling the model with the specified optimizer and loss function and metric functions
        model.compile(optimizer=adam, loss=loss_fit, metrics=[total_charge, AtomA_sum, AtomB_sum, q, rmse])

        # Defining an EarlyStopping callback
        early_stopping = EarlyStoppingWithBestWeights(monitor='loss', patience=stopping_patience, verbose=0, mode='min')

        # Training the model
        model_history = model.fit(coor_train, rho_train, batch_size = batch_s, epochs=epoch_num, verbose=0, callbacks=[early_stopping])

        # Collecting Training results:
        training_loss = model_history.history['loss']
        q_per_epoch = model_history.history['q']
        rmse_ep = model_history.history['rmse']

        

        # Create an array of epoch numbers
        epochs = range(1, len(training_loss) + 1)
        epoch_len = len(epochs)

        # Using the trained weights to get the AIM predicted density
        rho_pred = model.predict(coor_test, verbose=0)
        
        
        #Calculate the electric dipole moment at the origin by integrating
        #an electron density file with columns:  x  y  z  rho
        
        #    p_i = -e * sum[ r_i * rho(r) * dV ]
        
        #Coordinates in Bohr (Gaussian cubegen output), density in electrons/Bohr^3.
        #Output in Debye (1 e*Bohr = 2.54175 D).
        
        
        # Extract coordinates and density
        x = coor_test[:, 0]
        y = coor_test[:, 1]
        z = coor_test[:, 2] + shift
        

        
        # Voxel volume from grid spacing
        dx = np.unique(x)[1] - np.unique(x)[0]
        dy = np.unique(y)[1] - np.unique(y)[0]
        dz = np.unique(z)[1] - np.unique(z)[0]
        dV = dx * dy * dz
        
        # Dipole: p = -e * integral( r * rho dV )
        BOHR_TO_DEBYE = 2.54175
        w_train = rho_train * dV  # electrons per voxel
        px_train, py_train, pz_train = -(x * w_train).sum(), -(y * w_train).sum(), -(z * w_train).sum()
        
        
        w_pred = rho_pred * dV  # electrons per voxel
        px_pred, py_pred, pz_pred = -(x * w_pred).sum(), -(y * w_pred).sum(), -(z * w_pred).sum()
        
        
        # Electronic dipole (no sign yet)
        dipole_train = np.array([px_train, py_train, pz_train]) * BOHR_TO_DEBYE
        
        dipole_pred = np.array([px_pred, py_pred, pz_pred]) * BOHR_TO_DEBYE
        
        dipole_pred_lst.append(dipole_pred)
        dipole_train_lst.append(dipole_train)
        


        # Store trained weights:
        weights = model.layers[0].get_weights()
        flat_weights = [w[0] for w in weights]

        for i in range(N_A):
            wA_hist[i].append(flat_weights[i])

        for i in range(N_B):
            wB_hist[i].append(flat_weights[N_A + i])

        # Redefine initial weights so that next value of R starts with them
        ini_weights = np.array(flat_weights)


        # Calculate effective and total charges
        A_Eff = sum(
            flat_weights[i] * state["charge"]
            for i, state in enumerate(ATOM_A_STATES)
        )

        B_Eff = sum(
            flat_weights[N_A+i] * state["charge"]
            for i, state in enumerate(ATOM_B_STATES)
        )

        Tot_Ch = A_Eff + B_Eff

        
        # Print converged weights
        print("Given R:", R)
        print()

        for i, state in enumerate(ATOM_A_STATES):
            print(
                f"{AtomA}{state['name']} Weight:",
                flat_weights[i]
            )

        print("Sum:", sum(flat_weights[:N_A]))
        print()

        for i, state in enumerate(ATOM_B_STATES):
            print(
                f"{AtomB}{state['name']} Weight:",
                flat_weights[N_A+i]
            )

        print("Sum:", sum(flat_weights[N_A:]))

        print(f"\n{AtomA} Eff. Charge", A_Eff )
        print(f"{AtomB} Eff. Charge", B_Eff)
        print("\nTotal Charge: ", Tot_Ch)
        
        print(f"\nElectrons integrated (actual): {w_train.sum():.4f}")
        print(f"actual px = {px_train * BOHR_TO_DEBYE:+.4f} D")
        print(f"actual py = {py_train * BOHR_TO_DEBYE:+.4f} D")
        print(f"actual pz = {pz_train * BOHR_TO_DEBYE:+.4f} D")
        print(f"|p| = {np.linalg.norm([px_train, py_train, pz_train]) * BOHR_TO_DEBYE:.4f} D\n")
        
        print(f"\nElectrons integrated (pred): {w_pred.sum():.4f}")
        print(f"pred px = {px_pred * BOHR_TO_DEBYE:+.4f} D")
        print(f"pred py = {py_pred * BOHR_TO_DEBYE:+.4f} D")
        print(f"pred pz = {pz_pred * BOHR_TO_DEBYE:+.4f} D")
        print(f"|p| = {np.linalg.norm([px_pred, py_pred, pz_pred]) * BOHR_TO_DEBYE:.4f} D\n")
        
        
        # crude check on raw data (non-uniform grid, but diagnostic)
        dx = np.mean(np.diff(np.unique(x)))
        dy = np.mean(np.diff(np.unique(y)))
        dz = np.mean(np.diff(np.unique(z)))
        
        Ne_raw = np.sum(rho_train) * dx * dy * dz
        print(f"Electron count (raw approx): {Ne_raw}")


        rmse_val = np.sqrt( (1/len(rho_pred) * np.sum( (rho_pred-rho_test)**2 ) ))
        print("RMSE: ", rmse_val)
        
        try:
            for i in range(epoch_len):
                epoch_data.append([(i+1), training_loss[i], q_per_epoch[i], rmse_ep[i]])
        except Exception as e:
            print("Epoch data Error:")
            print(e)
        
        
        # Store results in mas_data
        row = [
            key_str,
            R,
            batch_s,
            lr,
            dcy
        ]

        row.extend(flat_weights)

        row.extend([
            A_Eff,
            B_Eff,
            Tot_Ch,
            rmse_val,
            epoch_len
        ])

        mas_data.append(row)

        df = pd.DataFrame(mas_data)
        df.to_csv(f'mas_data_T{trial}_{molecule}.csv')


        # Plot of loss as a function of epochs for this value of R
        plt.plot(epochs, training_loss, label='Training Loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.title('Epoch vs. Loss, {}A'.format(key_str))
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_EpochLoss_{key_str}A.png')
        plt.close()

        
        # Plot of charge on atom A as a function of epochs for this value of R
        plt.plot(epochs, q_per_epoch, label='Charge')
        plt.xlabel('Epochs')
        plt.ylabel(f'Charge on {AtomA}')
        plt.title(f'Epoch vs. {AtomA} Eff. Charge, {key_str}A')
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_EpochEffCharge_{key_str}A.png')
        plt.close()



        # Plot of RSME as a function of epochs for this value of R
        plt.plot(epochs, rmse_ep, label='RMSE')
        plt.xlabel('Epochs')
        plt.ylabel('RMSE of fit')
        plt.title(f'Epoch vs. RMSE, {key_str}A')
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_EpochRMSE_{key_str}A.png')
        plt.close()

        # Save training data for this value of R in a csv file
        df_sub = pd.DataFrame(epoch_data)
        df_sub.to_csv(f'T{trial}_{molecule}_epoch_data_{key_str}A.csv')
        

    except Exception as e:
        print("Error occured with distance {}A:".format(key_str))
        print(e)
        
        
    try:

        # ---------------------------------------
        # PLOTTING FINAL DATA AS A FUNCTION OF R
        # ---------------------------------------
        #     
        
        
        R_lst = []
        A_Eff_lst = []
        B_Eff_lst = []
        Ch_lst = []
        rmse_lst = []
        epoch_lst = []
        
        for line in mas_data[1:]:
            R_lst.append(line[1])
            A_Eff_lst.append(line[-5])
            B_Eff_lst.append(line[-4])
            Ch_lst.append(line[-3])
            rmse_lst.append(line[-2])
            epoch_lst.append(line[-1])
        
        size = 4
        max_R = (R_lst[0]) + 1
        min_R = round(R_lst[-1]) - 1
        
        
        
        
        # Plot of Effective Charge on both atoms and Total Charge
        plt.plot(R_lst, A_Eff_lst, '-o', markersize = size, label=f"{AtomA} Eff Ch.")
        plt.plot(R_lst, B_Eff_lst, '-o', markersize = size, label=f"{AtomB} Eff Ch.")
        plt.plot(R_lst, Ch_lst, '-o', markersize = size,  label="Total Charge")
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Charge")
        plt.title(f'{level_of_theory} {molecule} Charge vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.yticks(np.arange(-1,1.2,0.05), minor=True)
        plt.yticks(np.arange(-1,1.2,0.1))
        
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_SepComplCharge.png')
        plt.close()
        
        
        
        
        #Plot of Effective charge only for Atom A
        plt.plot(R_lst, A_Eff_lst, '-o', markersize = size, label=f"{AtomA} Eff Ch.")
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Charge")
        plt.title(f'{level_of_theory} {AtomA} Charge vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.yticks(np.arange(0,1.2,0.05), minor=True)
        plt.yticks(np.arange(0,1.2,0.1))
        
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_SepCharge.png')
        plt.close()
        
        
        
    
        
        # Plot of number of epochs required for each value of R
        plt.plot(R_lst, epoch_lst, '-o', markersize = size, label="epochs")
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Number of Epochs")
        plt.title(f'{level_of_theory} {molecule} Epoch Number vs. Separation, Lr={lr}, Dcy={dcy}')
        
        plt.legend()
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.savefig(f'T{trial}_{molecule}_SepEpochLen.png')
        plt.close()
        
        
        
        
        
        
        # Plot of RSME for each value of R
        plt.plot(R_lst, rmse_lst, '-o', markersize = size, label="epochs")
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("RMSE")
        plt.title(f'{level_of_theory} {molecule} RMSE vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_SepRMSE.png')
        plt.close()
        
        
        
    
        
        #Plot of Weights for Atom A
        for i, state in enumerate(ATOM_A_STATES):
            plt.plot(
                R_lst,
                wA_hist[i],
                '-o',
                markersize=size,
                label=f"{AtomA}{state['name']}"
            )
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Weight")
        plt.title(f'{level_of_theory} {AtomA} Weights vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.yticks(np.arange(0,1,0.05), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.yticks(np.arange(0,1.1,0.1))
        
        plt.legend()
        plt.savefig(f'T{trial}_{AtomA}_SepWeights.png')
        plt.close()
        
        
        #Plot of Weights for Atom B
        
        for i, state in enumerate(ATOM_B_STATES):
            plt.plot(
                R_lst,
                wB_hist[i],
                '-o',
                markersize=size,
                label=f"{AtomB}{state['name']}"
            )
        
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Weight")
        plt.title(f'{level_of_theory} {AtomB} Weights vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.yticks(np.arange(0,1,0.05), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.yticks(np.arange(0,1.1,0.1))
        
        plt.legend()
        plt.savefig(f'T{trial}_{AtomB}_SepWeights.png')
        plt.close()
        
        
        
        
        
        
        #Plot of All weights
        
        for i, state in enumerate(ATOM_A_STATES):
            plt.plot(
                R_lst,
                wA_hist[i],
                '-o',
                markersize=size,
                label=f"{AtomA}{state['name']}"
            )
        
        for i, state in enumerate(ATOM_B_STATES):
            plt.plot(
                R_lst,
                wB_hist[i],
                '-o',
                markersize=size,
                label=f"{AtomB}{state['name']}"
            )
        

        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Weight")
        plt.title(f'{molecule} Weights vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R,max_R,0.5), minor=True)
        plt.yticks(np.arange(0,1,0.05), minor=True)
        plt.xticks(np.arange(min_R,max_R,1.0))
        plt.yticks(np.arange(0,1.1,0.1))
        
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_SepWeights.png')
        plt.close()
        
        
        
        
        # Dipole Estimate Calculation for as a function of R:
        
        
        dipole_pred_arr  = np.array(dipole_pred_lst)   # shape (n_R, 3)
        dipole_train_arr = np.array(dipole_train_lst)
        
        #dipole_pred_arr[:, 2] *= -1  # flip z-component
        #dipole_train_arr[:, 2] *= -1  # flip z-component
        
        for i, component in enumerate(['x', 'y', 'z']):
            plt.plot(R_lst, dipole_pred_arr[:, i],  '-o', markersize=size, label=fr"Predicted $\mu$_{component}")
            plt.plot(R_lst, dipole_train_arr[:, i], '-o', markersize=size, label=fr"Training $\mu$_{component}")
        
            plt.xlabel("Separation (a.u.)")
            plt.ylabel("Dipole Moment")
            plt.title(fr'{level_of_theory} {molecule} Dipole Moment $\mu$_{component} vs. Separation, Lr={lr}, Dcy={dcy}')
            plt.xticks(np.arange(min_R, max_R, 0.5), minor=True)
            plt.xticks(np.arange(min_R, max_R, 1.0))
            plt.legend()
            plt.savefig(f'T{trial}_{molecule}_SepDipole_{component}.png')
            plt.close()
            
            
            
            
        # directory name like "0.90A"
        r_dir_re = re.compile(r'^(\d+(\.\d+)?)A$')
        
        def read_file_text(fname):
            with open(fname, "r", errors="replace") as f:
                return f.read()
        
        def extract_dipole(text):
            dipz_re = re.compile(
                r'Dipole moment \(field-independent basis, Debye\):\s*.*?Z=\s*([-+]?\d+\.\d+)',
                re.DOTALL | re.IGNORECASE
            )
            dips = dipz_re.findall(text)
            return float(dips[-1]) if dips else None
        
        
        system = molecule
        
        
        plt.figure()
        
        direc = f"/easley/scratch/users/sol-sam/FineMolecularDensities/{molecule}/{level_of_theory}/Fine/{molecule}_{level_of_theory}_{basis_set}_seq/"
        level = f"{level_of_theory}"
        data = []  # (R, dipole)
    
        for name in sorted(os.listdir(direc)):
            dpath = os.path.join(direc, name)
            if not os.path.isdir(dpath):
                continue
    
            m = r_dir_re.match(name)
            if not m:
                continue
    
            R = float(m.group(1))
            expected_log = os.path.join(dpath, f"{system}_{name}.log")
    
            log_file = expected_log if os.path.isfile(expected_log) else None
    
            if log_file is None:
                for fn in os.listdir(dpath):
                    if fn.endswith(".log"):
                        log_file = os.path.join(dpath, fn)
                        break
    
            if log_file is None:
                continue
    
            text = read_file_text(log_file)
            dip = extract_dipole(text)
    
            if dip is not None:
                data.append((R, dip))
    
        data.sort(key=lambda x: x[0])
    
        Rs = [x[0] * 1.8897259886 for x in data]
        dips = [x[1] for x in data]  # flip sign
    
        plt.plot(Rs, dips, marker='o', markersize=size, label=f"{level} Gaussian Value")
        
        
        plt.axvline(x=1.128 * 1.8897259886, linestyle='--', label='Equilibrium separation')
        
        
        plt.plot(R_lst, dipole_pred_arr[:, 2],  '-o', markersize=size, label=fr"Predicted $\mu$_z")
        plt.plot(R_lst, dipole_train_arr[:, 2], '-o', markersize=size, label=fr"Training $\mu$_z")
    
        plt.xlabel("Separation (a.u.)")
        plt.ylabel("Dipole Moment")
        plt.title(fr'{level_of_theory} {molecule} Dipole Moment $\mu$_z vs. Separation, Lr={lr}, Dcy={dcy}')
        plt.xticks(np.arange(min_R, max_R, 0.5), minor=True)
        plt.xticks(np.arange(min_R, max_R, 1.0))
        plt.xlim(min_R, max_R)
        plt.legend()
        plt.savefig(f'T{trial}_{molecule}_SepDipole_GaussianCompare.png')
        plt.close()
        
    except Exception as e:
        print("Plotting Error occured with distance {}A:".format(key_str))
        print(e)
        
        
    sep_end = time.perf_counter()
    sep_elapsed = sep_end - sep_start
    
    minutes = int(sep_elapsed // 60)
    seconds = int(sep_elapsed % 60)
    
    print(f"Elapsed time for {key_str} A: {minutes} min {seconds} s")



total_end = time.perf_counter()
total_elapsed = total_end - total_start
total_minutes = int(total_elapsed // 60)
total_seconds = int(total_elapsed % 60)

print(f"Total runtime: {total_minutes} min {total_seconds} s")

# Job Complete Output
print("\n--------------------------------------------------------------------")
print("Job Complete!\n")
current_datetime = datetime.now()
today = date.today().strftime("%m/%d/%y")
cur_time = current_datetime.time()
print(f"Sol S., {today}  {cur_time.hour}:{cur_time.minute}:{cur_time.second}\n")


