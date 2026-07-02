*******************************************
README File for Essential Support Files For MANN implementation 
Last Updated: 7/2/2026
*******************************************

Scripts/Files within this directory are essential support files for the 
implementation of the RBF_NN.py script. This means that these scripts are 
called by the main script (RBF_NN.py) and are used "under the hood" for 
processing the radial basis function distributions for given atomic states
needed by the RBF-NN. 

These files MUST be placed in the same directory as RBF_NN.py when the main
script is started. These files do not need to be edited between changes in 
the molecular system being studied or changes in the number of states being 
requested. These files would need to be edited if a new state is added to the 
library of states callable by the RBF-NN or if a new RBF Model is added 
(currently M5, M6, and M7) are callable (see below for further description).



Files within this directory include:

* RBF_ModelParameters.xlsx
	- excel spreadsheet containing model parameters for all states fitted 
	  by the MATLAB RBF model fitting routine 
	  (see \RadialBasisFunctions_MatLabFitting directory for more details 
	  on process).
        - See work in [1] and [2] for details on the theory and naming 
	  convention behind the parameters.
	- These parameter values are then pulled by params_loader.py for use in
	  the main RBF-NN script.
	- In the case that a new atomic state is fitted and is to be studied with
	  the RBF-NN script, the fitted parameters of the new state will need to 
	  be added as a new row in the excel spreadsheet following the same format
	  as currently implemented states. 

* params_loader.py
	- This python script is called by the main RBF_NN.py as a library script
	  in order to pull fitted RBF model parameters from RBF_ModelParameters.xlsx
	- The collected model parameters are then used by the rbf_functions.py script
	  in order to determine the full radial density function for a state. 
	- In the case that a new atomic state is to be added to the RBF-NN script,
	  the directory "_STATE_ALIAS" must be edited to add the name of the new 
	  state following the format of previous states.

* rbf_functions.py
	- This python script is called by the main RBF_NN.py as a library script
	  in order to output the radial density distribution for a specific state.
	- The script contains functions for Models M5, M6, and M7 as according to
	  work described in [1] and [2]. The script also contains the function
	  "hydrogen_density(self, coor, R, n, l)" to output radial densities for
	  the ground and excited states of hydrogen given the quantum numbers (n,l)
	- See work in [1] and [2] for details on the theory behind these models
	  (M5, M6, M7) and their usage.
	- If adding a new model to be used by the RBFNN, a new function of similar
	  format to M5, M6, M7 would need to be added to this script that would be 
	  imported accordingly to the RBF_NN.py script.

References: 
[1] Godwin Amo-Kwao. Radial basis densities and the density functional-based 
    atom-in molecule: Designing charge-transfer potentials. Ph.D. dissertation, 
    University of New Mexico, 2020
[2] Godwin Amo-Kwao, Sol Samuels, and Susan R. Atlas. Radial basis function 
    electron densities with asymptotic constraints, manuscript 2026