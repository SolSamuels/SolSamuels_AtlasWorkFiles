
######################################
# rbf_functions.py by Sol Samuels
# last edited, 7/2/2026
######################################

# Script to implement the radial basis function models (M5, M6, M7) for 
# a given state in the RBF-NN keras layer. Function for the ground and excited states
# of neutral hydrogen is also included, hydrogen_density(self, coor, R, n, l), which determines
# the radial density distribution given the quantum numbers (n,l) for the state.
# This is a support script for RBF_NN.py (head script for RBF-NN implementation)
# No edits should be needed here to run, even if switching directory location / user. However, this 
# script must be in the same directory as the head script RBF_NN.py to run.


import tensorflow as tf
import matplotlib.pyplot as plt
from scipy.integrate import trapz


# Model M5, used for neutral and cation states
def M5(self, coor, R, params):
    
    r = tf.sqrt((coor[:,0])**2 + (coor[:,1])**2 + (coor[:,2] - R)**2)
    
    A_term = params['A0'] * tf.exp(-2 * params['Z'] * tf.abs(r))
    B_term = params['B0'] * tf.abs(r)**(2*params['beta']) * tf.exp(-2 * params['alpha'] * tf.abs(r))
    C_term = params['C0'] * (params['mm'] - params['Z'] * tf.abs(r))**2 * tf.exp(-2 * params['gamma'] * tf.abs(r))
    D_term = params['D0'] * (params['Z'])**2 * (tf.abs(r))**2 * tf.exp(-2 * params['eta'] * tf.abs(r))
    M5 = A_term + B_term + C_term + D_term
    return M5


# Model M6, used for anion states
def M6(self, coor, R, params): 
    
    r = tf.sqrt((coor[:,0])**2 + (coor[:,1])**2 + (coor[:,2] - R)**2)
    
    beta_eff = params['beta'] + params['bb'] * tf.abs(r)

    A_term = params['A0'] * tf.exp(-2 * params['Z'] * tf.abs(r))
    B_term = params['B0'] * tf.abs(r)**(2*beta_eff) * tf.exp(-2 * params['alpha'] * tf.abs(r))
    C_term = params['C0'] * (params['mm'] - params['Z'] * tf.abs(r))**2 * tf.exp(-2 * params['gamma'] * tf.abs(r))
    D_term = params['D0'] * (params['Z'])**2 * (tf.abs(r))**2 * tf.exp(-2 * params['eta'] * tf.abs(r))
    M6 = A_term + B_term + C_term + D_term
    return M6


# Model M7, used for excited states
def M7(self, coor, R, params): 
    
    r = tf.sqrt((coor[:,0])**2 + (coor[:,1])**2 + (coor[:,2] - R)**2)
    
    beta_eff = params['beta'] + params['bb'] * tf.abs(r)

    A_term = params['A0'] * tf.exp(-2 * params['Z'] * tf.abs(r))
    B_term = params['B0'] * tf.abs(r)**(2*beta_eff) * tf.exp(-2 * params['alpha'] * tf.abs(r))
    C_term = params['C0'] * (params['mm'] - params['Z'] * tf.abs(r))**2 * tf.exp(-2 * params['gamma'] * tf.abs(r))
    D_term = params['D0'] * (params['Z'])**2 * (tf.abs(r))**2 * tf.exp(-2 * params['eta'] * tf.abs(r))
    E_term = params['E0'] * (params['pp'] - (params['qq']*params['Z']*tf.abs(r)) + 
              (params['tt']*params['Z']*tf.abs(r))**2)**2 * tf.exp(-2 * params['delta'] * tf.abs(r))
    M7 = A_term + B_term + C_term + D_term + E_term
    return M7


# function to determine density distributions for states of neutral hydrogen
# This function was originally coded in MATLAB by Susan R. Atlas and Steven Boyd and
# translated to python by Sol Samuels

def hydrogen_density(self, coor, R, n, l):
    r = tf.sqrt((coor[:,0])**2 + (coor[:,1])**2 + (coor[:,2] - R)**2)
    
    if n < 1:
        raise ValueError('Error: invalid value of n; Exiting...')
    elif n > 4:
        raise ValueError('Value of n exceeds supported values; Exiting...')
    
    Z = 1.0  # Nuclear charge (1 for H, but keep in code for generalizability to other hydrogenic wavefunctions.
    rscal = 2 * Z / n  # Scale factor for distance (see: http://en.citizendium.org/wiki/Hydrogen-like_atom, referencing Pauling and Wilson)
    Rnorm = Z**1.5  # Normalization factor for Rnl

    # Scale the radius
    rmatscal = r * rscal # scaled r values for specified value of n (denoted by rho_n in most texts).
    expfac = tf.exp(-rmatscal * 0.5) # exponential prefactor common to all Rnl

    # compute radial wavefunction x angle-averaged angular wfn (Anglenorm) for quantum numbers (nl) and store in fMat
    if n == 1 and l == 0:  # R10
        fmat = Rnorm * expfac * 2.0
    elif n == 2 and l == 0:  # R20
        fmat = Rnorm * expfac * (0.5 / tf.sqrt(2.0)) * (2.0 - rmatscal)
    elif n == 2 and l == 1:  # R21
        fmat = Rnorm * expfac * (0.5 / tf.sqrt(6.0)) * rmatscal
    elif n == 3 and l == 0:  # R30
        fmat = Rnorm * expfac * (1.0 / (9.0 * tf.sqrt(3.0))) * (6.0 - 6.0 * rmatscal + rmatscal**2)
    elif n == 3 and l == 1:  # R31
        fmat = Rnorm * expfac * (1.0 / (9.0 * tf.sqrt(6.0))) * (4.0 - rmatscal) * rmatscal
    elif n == 3 and l == 2:  # R32
        fmat = Rnorm * expfac * (1.0 / (9.0 * tf.sqrt(30.0))) * rmatscal**2
    elif n == 4 and l == 0:  # R40
        rmatscalsq = rmatscal**2
        fmat = Rnorm * expfac * (1.0 / 96.0) * (24.0 - 36.0 * rmatscal + 12.0 * rmatscalsq - rmatscalsq * rmatscal)
    elif n == 4 and l == 1:  # R41
        rmatscalsq = rmatscal**2
        fmat = Rnorm * expfac * (1.0 / (32.0 * tf.sqrt(15.0))) * (20.0 - 10.0 * rmatscal + rmatscalsq) * rmatscal
    elif n == 4 and l == 2:  # R42
        rmatscalsq = rmatscal**2
        fmat = Rnorm * expfac * (1.0 / (96.0 * tf.sqrt(5.0))) * (6.0 - rmatscal) * rmatscalsq
    elif n == 4 and l == 3:  # R43
        rmatscalcub = rmatscal**3
        fmat = Rnorm * expfac * (1.0 / (96.0 * tf.sqrt(35.0))) * rmatscalcub
    else:
        raise ValueError(f'Unsupported quantum numbers n={n}, l={l}')
        
    pi = 3.141592653589793

    # rhoMat (total density) is the radial wfn squared, times 1/(4*pi) for the angle-averaged (theta,phi) wfn-sqd component.
    rho = (fmat**2) / (4.0 * pi)

    return rho


# In[ ]:






