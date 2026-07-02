%% Script to call densityFitting.m 
% ========================================================================
% This script sets up and calls the densityFitting.m 

clc;
clear;
close all;

%% ------------------- USER INPUT SECTION -------------------
% Atomic state: (should match filenames)
atomicstate = 'O_S4';


% Input filenames (exact same variable names as original code expects)
input_param  = sprintf('SolGen_Collected_Input_Files/%s_param.dat', atomicstate);     % Parameter file specific to atom/anion
input_bcwgt  = sprintf('SolGen_Collected_Input_Files/%s_input_bcwgt.dat', atomicstate);   % Boundaries and weights file 
% g09_data     = sprintf('Exct_Collected_Input_Files/%s_g09_zcut.dat', atomicstate);  % Gaussian zcut data file
dBardata     = sprintf('SolGen_Collected_Input_Files/%s_dBar.dat', atomicstate);      % Sphericalized dBar data file

% Data flag: either 'g09' or 'dBar' depending on the data type to be fitted
dataFlag     = 'dBar';  

%% ------------------- CHECK FILE EXISTENCE -------------------

disp('------------------------------------------------------');
disp('Checking required input files for density fitting...');
if ~isfile(input_param)
    error('Error: Input parameter file %s not found.', input_param);
end
if ~isfile(input_bcwgt)
    error('Error: Input weights/bounds file %s not found.', input_bcwgt);
end
if ~isfile(dBardata)
    error('Error: dBar radial data file %s not found.', dBardata);
end
disp('All input files located successfully.');
disp('------------------------------------------------------');

%% ------------------- CALL FITTING ROUTINE -------------------

disp('Starting density fitting...');
try
    % Call the function exactly as expected
    densityFitting_m7(input_param, input_bcwgt, dBardata, dataFlag);
    disp('Density fitting completed successfully.');
catch ME
    % Catch and display error if fitting fails
    disp('An error occurred during density fitting:');
    disp(getReport(ME));
end

%% ------------------- FINISHED -------------------
disp('Script execution completed.');
disp('------------------------------------------------------');
