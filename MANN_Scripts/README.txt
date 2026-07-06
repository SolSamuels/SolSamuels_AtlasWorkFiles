*******************************************
README File for MANN Work by Sol Samuels Completed
within the Atlas research group, Summer 2023 - Summer 2026

Last Updated: 7/6/2026
*******************************************

This repository contains the most essential files/computational programs
needed to continue or recreate work done for the Radial Basis Function Neural Network
(RBFNN) project (also referred to as "Molecules as Neural Networks" (MANN)).


Files/Subdirectories within this directory include:

* \RBFNN :
	- This subdirectory contains the main script for the RBF-NN program, 
          titled "RBF_NN.py", as well as an example slurm submission script 
	  "RBF_NN.sh".
	- For more information for these files please see the subdirectory 
	  readme file \RBFNN\README.txt

* \EssentialSupportFiles :
	- This subdirectory contains the essential support scripts/files needed
	  to run RBF_NN.py (found in \RBFNN). These support files need to be within
	  the same directory as RBFNN.py in order for that script to run.
	- For more information for these files please see the subdirectory readme file
	  \EssentialSupportFiles\README.txt

* Environment_Setup.txt
	- This file contains ESSENTIAL instructions for how to setup an environment
          in CARC such that all scripts in this directory can run. This includes how
	  to include tensorflow, scipy, and other libraries in your environment.
        - This should be your next step.
	
* \ExampleRun :
	- This subdirectory contains an example setup and example output for running
	  RBF_NN.py as well as instructions for how to run.

