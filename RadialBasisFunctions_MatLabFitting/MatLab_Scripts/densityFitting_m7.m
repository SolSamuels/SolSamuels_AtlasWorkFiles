% Density fitting code for neutral and charged atomic states.
%
% Author: Godwin Amo-Kwao (initial version, 2018); unified version for neutral, cation, and anion states, with
% added documention (2019; 2020).
%
% Modified by S. R. Atlas to add options for different dBar fitting strategies.
%
% Version 1.1, GAK, 7/12/19.      Initial grand-unified version.
% Version 1.2, SRA, 7/14/19.    * Fixed Angstrom to Bohr conversion and made it a defined constant.
%                               * Clarified some of the documentation and formatting, and fixed typos
%                                 on output/plots.
%                               * Fixed sanity check on \rho(0) to express in terms of relative, 
%                                 not absolute, tolerance.
%                               * Fixed computation of B0_dBar to use dBar values and radial grid
%                                 instead of G09 data and zcut grid.
%
% Version 1.3, GAK/SRA 7/17/19. * Fixed additional typos, and added missing 'hold on' for dBar plot
%                               * Fixed factors of 4*pi and defined 4*pi as a constant
%                               * Got rid of B0_dBar and B0_G09--- not needed since fitting only to dBar data.
%                                 Report A0_G09 and A0_dBar in output, but set A0 = A0_G09 for fit,
%                                 to satisfy Kato constraint at nucleus.
%
% Version 1.4, SRA 7/19/19.     * Renamed nbo_param to more generic input_param file
%                               * Fixed output column labels to be more informative.
%                               * Added new parameter, 'model' (hardwired for now, but to be read in
%                                 via input_param file), in preparation for testing Model 2 (beta as fitting parameter rather than
%                                 input), and Model 3 (Wang & Parr), in addition to analytic Model 1.
%                               * Put in placeholders for Model = 2 and Model = 3 fitting options for electron density.
%                               * Changed output density plots to be the more useful ln \rho rather than \rho.
%                               * Renamed den_analyt_int etc. to model_fit to be generic for all 3 models
%                               * Changed hardwired plotting range for ln densities to xlim [0 3]
%
% Version 1.5, GAK/SRA 7/26/19. * Added additional log file at top of code, to enable formatted printing of key parameters
%                               * Added initialization of B0 and alpha using Model 1 as input to Model 2
%                                 fmincon optimization.
%                               * Implemented Model 2 as fmincon optimization of B0 and alpha parameters,
%                                 with constrained ranges as described in dissertation text,
%                                 and overall constraint for model to integrate to total # electrons (Ne)
%                                 incorporated within objective function.
%                               * Tried implementing weights within objective function to control relative importance of fit to density and integration of
%                                 model density to Ne (the weights sum to 1, currently hardwired within objective-function function.
%
% Version 1.6, SRA/GAK 8/26/19. * Implemented Model 3, with additional term to address 2s shell structure (in the spirit
%                                 of Karplus & Porter hydrogenic function.
%
% Version 1.7, SRA/GAK 9/2/19.  * Fit with RDF data in objective function.
%                               * Added weights to regions R1 - R4 in objective function.
%                          
% Version 1.8, SRA/GAK 9/11/19. * Added A0 as optimization parameter.
%                               * Added beta as optimization paremeter.
%                               * Added omputation of exact Kato condition Kato = (-\rho'(0) /(2*\rho(0))
%                                 as part of the optimization process, with corresponding contribution
%                                 to the objective function (Kato-1.0)^2, and weight w_5.
%                               * Added printout of comma-delimited .txt file to facilitate import into
%                                 run/input/output parameter Excel spreadsheet.
%                               * Added check that weights sum to 1.
%
% Version 1.9, SRA/GAK 9/15/19. * [Still to be added]. Handler to smoothly transition from 2p orbital density (long-range constraint)
%                                 appropriate for N=Z, to 2p orbital density (long-range constraint with N=Z 2p orbital
%                                 density, and extra electron far away), for N=Z+1.  This may only be needed only
%                                 for Z > 4, when 2p orbital is occupied in the neutral atom. Analogous handling
%                                 will be needed for transitioning the 2s orbital when Z <=4 and N = Z+1.
%                               * Documented code to allow easy toggling between fitting to G09 (rGrid, rho_g09) or
%                                 dBar (rdBar,rhodBar) data, for testing.
%                               * Fixed bug in computation of Kato (had "alpha" instead of "2Z" in both numerator
%                                 and denominator.
%
% Version 1.95, GAK 6/12/20.    * Added version check for readtable formatting of paramfile
%                               * Added data flag to calling sequence for selecting g09 or dBar
%                               % Made rGrid generic to apply to either g09 or dBar fitting 
%                               * Added some missing boundary conditions to file generated for Excel.
%                               * Code clean up.
%
% Version 1.96, GAK 7/09/20     * Added Model 4 that uses Zeff as a fitting parameter.
%                               * Added padcat routine to automate padding of Matlab 'vectors' of different lengths 
%                                 for plot compatibility (Model 2, 3, 4 comparison) by padding with NaNs                            
%
% Version 1.97,  GAK 7/18/20    * Added savePlots for saving high resolution cropped png/pdf images in atomic state directories.
% Version 1.98,b GAK 7/27/20    * Fixed w2/w3 to be consistent with F2/F3 in Model 3 & 4 objective functions.
%                               * w1, w2, w3, w4, w5 in this version corresponds to old w1, w3, w2, w4, w5.
%                                 This version fixes a bug that has existed at least since version 1.95 (and most
%                                 likely since version 1.7, when different weights for different regions were implemented.
%                               * Tested that Model 4 reduces to Model 3 when Zeff is fixed at Model 4 fitted value.
%
% Version 1.99, SRA 9/10/20     * Implemented Model 5 to accommodate anion fitting
%                               * Moved 'smallnum' definition earlier in
%                                 code, to appear prior to bifurcation to individual
%                                 models.
%                               * Code cleanup.
%
% Version 2.0, SRA 10/8/20      * Added Model 5 and corresponding Model5Objfn.m and density_model_5_eval_fn.m
%                                 Model 5 incorporates an additional term (with coefficient D0 and new exponential decay
%                                 constant eta) to add an H-like 2p radial density to the fitting function. 
%
%                                 Note: Model 5 has 10 fitting parameters, although with the extra 2p flexibility,
%                                 it may make sense to fix alpha, beta, and Zeff at their nominal values and only
%                                 vary 7 parameters when fitting with this model.
%
%                               * Added input parameters for Model 5 (D0, eta, and their upper and lower bounds) to input_bcwgt file, 
%                                 and for reading in here (see Model 5 section).
%                               * Fixed bug in computation of Kato (had "gamma" instead of "2*gamma" in both numerator
%                                 and denominator.) This affected Model 3, 4, and 5 objective function routines, which are 
%                                 corrected in this version.
%                               * Code cleanup
%
% Version 2.02, SRA 10/15/20    * Added Model 6 and corresponding objective function and evaluation function.  This adds one new parameter, bb, to 
%                               *    in order to model beta with a linear fit:  \beta = \beta_0 + bb* r, where 
%                               *    bb is a new input parameter, and \beta_0 is the standard \beta parameter of
%                               *    earlier models.  Set bb = 0.0 to default to previous models.   Model 6 is intended
%                               *    to help improve the fit at long range for anions, and interpolate between atomic and
%                               *    anionic long-range limits.
%                               * Changed initialization of beta to be betaLB instead of nominal -1.  
%                               *     betaLB can be set = -1 or the nominal beta from the neutral. 
%                               * Fixed bug in test for weights not = 1 (added absolute value). Code will now throw error and exit.
%                               * Modified computation of long-range function on Model6ObjFn for consistent handling of r=0 as
%                               *     elsewhere in code, and .* for bb multiplication (no impact on results). 
%                               * Added printout of effective beta at short- and long-range; also plot of (linear) beta as
%                               *      a function of r.
%
% Function arguments include the following three files:
%
% "input_param" file contains the parameters used for fitting.
%  Line 1: atomic species (string)
%  Line 2: highest occupied orbital energy value (eV) obtained from NBO analysis applied to G'09-computed results.
%          Alternatively, use experimentally-determined electron affinity.
%  Line 3: Atomic number; number of electrons
%  Line 4: Model selection (1, 2, 3, 4, 5, or 6).
%
% "g09_data" file contains zcut of 3D electron density as (z, \rho), generated from Gaussian '09's cubegen subroutine
%
% "dBardata" file contains radial density (r, \rho) generated by applying sphericalization code IntegrateDensity_fn to
%  density data from Gaussian '09's cubegen subroutine.
%
% "input_bcwgt" file contains information for fitting: initial guesses for parameters C0, mm, and gamma;
% lower and upper range constraints for fitting parameters (A0, B0, C0, mm, gamma, and Zeff); objective function 
% weights w1-w5; range parameters r1, r2 and cutoff rCut. 
% 
% Set dataFlag to 'g09' to fit to Gaussian '09 data or 'dBar' to fit to dBar data.
%
% Example invocation of the code:
% densityFitting('li0_param.dat','input_bcwgt.dat','li0_g09_zcut.dat','li0_dBar.dat','g09')

function densityFitting_m7(input_param,input_bcwgt,dBardata,dataFlag)
close all;
format long

% Toggle global search
useGlobalSearch = true;   % true = particleswarm + fmincon
                          % false = fmincon only

% Constants
AngstromToBohr = 1.0/0.529177210903; % NIST CODATA conversion,
                                     % https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
fourPi = 4.0*pi;

% set timestamp
dateStrTime = datestr(now,'yymmdd_HHMMSSFFF');
dateNumTime = datetime('now','Format','yyyyMMdd.HHmmssSSS');
dateNumTime=datenum(dateNumTime);

%% Read input_param file
% check matlab version v = ver('MATLAB'); This is for changes in readtable in ver 9.8 (not back-compatible)
if verLessThan('matlab','9.8')
    % Import first line of file for atomic state under consideration
    atomicstate = readtable(input_param, 'ReadVariableNames',false);
else
    % Import first line of file for atomic state under consideration
     atomicstate = readtable(input_param, 'ReadVariableNames',false, 'Format','auto');
end

% Convert to character type for naming output files
atomicstate = char(table2array(atomicstate(1,2)));

% Import NBO (natural bond orbital) parameters
input_param = readtable(input_param, 'HeaderLines',1);

% Convert to output type
input_param = table2array(input_param(:,2));

% Allocate nbo data array to column variable names
% Note: The electron affinity (EA) of a neutral atom is the same as the ionization
% potential (IP) of the corresponding negative ion.  The highest occupied, least
% negative natural bond orbital (NBO) single-particle energy of the negative ion
% approximates the IP of the negative ion.

nbo_enrgy = input_param(1,1);    % Either experimental or from ab initio calculation
Z     = input_param(2,1);        % Z = atomic_number
inZ   = Z;                       % Used for documenting in output to Excel 
Ne    = input_param(3,1);        % Number of electrons
model = input_param(4,1);         % Model selection (1, 2, 3, 4, or 5)
eps   = abs(nbo_enrgy);          % epsilon (highest occupied NBO energy)
alpha = sqrt(2*eps);             % 1st long-range density parameter
beta  = ((Z-Ne+1)/alpha - 1);    % 2nd long-range density parameter (may be replaced by optimized fitting parameter depending on model chosen)
%% End input data read

modelNam = num2str(model);

% Open log file for writing; name with atomic state and time
% Specify output folder
plot_path = char(sprintf('./outputs/%s/',atomicstate));
% check if directory exists; if not, create
if ~exist(plot_path, 'dir')
    mkdir(plot_path)
end

outLog1 = char(sprintf('logfile_m%s_%s.txt', modelNam, dateStrTime));
outLog1 = strcat(plot_path,outLog1);
fid1    = fopen(outLog1,'wt');
outLog2 = char(sprintf('table_m%s_%s.txt', modelNam, dateStrTime));
outLog2 = strcat(plot_path,outLog2);
fid2    = fopen(outLog2,'wt');

% test to make sure file opened successfully
assert(fid2>0)

outLog3 = char(sprintf('excelData_m%s_%s.txt', modelNam, dateStrTime));
outLog3 = strcat(plot_path,outLog3);
outLog4 = char(sprintf('plotData_m%s_%s.txt', modelNam, dateStrTime));
outLog4 = strcat(plot_path,outLog4);



% Read in dBar data.  Note that dBar includes factor of 4*pi.
% rdBar is the array of r values output by the sphericalization code

dBardata = dlmread(dBardata);
rdBar    = dBardata(:,1);
rhodBar  = dBardata(:,2);

% Divide by 4*pi to get spherically-averaged density at nucleus comparable to G'09 total density there.
% Due to numerical errors in angular integration near r=0, set A0 for dBar fitting to be exactly equal to \rho(0) from G09, i.e. A0_G09, for testing.
A0_dBar  = dBardata(1,2)/(4*pi);   % 1st row, 2nd column for rho_max_dBar

% Set variables depending on whether fitting g09 or dBar data
if strcmp(dataFlag,'g09')
    fprintf('Warning! This version of code cannot fit zcut data');
    
elseif strcmp(dataFlag,'dBar')
%-- fit to dBar data --
    A0     = A0_dBar;
    inR    = rdBar;
    inRho  = (rhodBar/fourPi);
%-- only used for models 3, 4, 5, 6 --
    ini_A0 = A0_dBar;
else
    errMsg1 = 'ERROR: dataFlag must be set to dBar';
    error(errMsg1)
end

if(model == 1)
%% ------------------------------ Model 1 -------------------------------------
% Evaluate B_0 using analytic model constraints at nucleus (cusp condition) and long range
% (see GAK dissertation text for details).

B0 = intgrl_B0_num(A0, Z, alpha, Ne);

%% ------------------------------ Model 2 -------------------------------------
% Numerical fit to determine B0 and beta (see GAK dissertation text for details)

elseif(model == 2)

% Used in Model 2 (B0LB, alphaLB) only
smallNum = 1.e-6;
fprintf(fid1,'Small number = %f.\n\n',1.e-6);
    
% Define objective function for fminsearch as a function of x
fun2 = @(x)Model2ObjFn(x, inR, inRho, A0, Z, Ne);

% Run analytic Model 1 to compute initial guess for B0 parameter.
B0  = intgrl_B0_num(A0, Z, alpha, Ne);

fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha from Model 1 = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);

% initial parameter guess for fminsearch
x0  = [B0,alpha];
A   = [];
b   = [];
Aeq = [];
beq = [];

B0LB = smallNum;
B0UB = 0.05;

alphaLB = smallNum;
alphaUB = 0.64;

fprintf(fid1,'Upper and lower bound constraints for B0:    [%f,%f]\n',B0LB, B0UB);
fprintf(fid1,'Upper and lower bound constraints for alpha: [%f,%f]\n',alphaLB, alphaUB);

lb = [B0LB, alphaLB];
ub = [B0UB, alphaUB];

% From previous version using fminsearch instead of fmincon --------------------
% Optimize B0 and alpha parameters using fminsearch
% options  = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-6);
% bestx    = fminsearch(fun,x0,options);
% ------------------------------------------------------------------------------
bestx      = fmincon(fun2,x0,A,b,Aeq,beq,lb,ub);

B0         = bestx(1);
alpha      = bestx(2);
beta       = ((Z-Ne+1)/alpha - 1);

fprintf(fid1,'Optimized B0 from Model 2 = %f\n',B0);
fprintf(fid1,'Optimized alpha from Model 2 = %f\n',alpha);
fprintf(fid1,'Optimized beta computed from Model 2 alpha = %f\n',beta);

assert(fclose(fid1)==0);               % check for successful file close

%% ------------------------------ Model 3 -------------------------------------
% Numerical fit with addition of hydrogen-like 2s orbital for Period 2 elements
% (see GAK dissertation text for details)

elseif(model == 3)

% Import boundary conditions and weights
input_bcwgt = readtable(input_bcwgt, 'HeaderLines',0);

% Convert to desired output format (preserves strings as strings and
% numbers as numbers)
input_bcwgt = table2array(input_bcwgt(:,2));

% Run analytic Model 1 to compute initial guess for B0.  Initial guesses for
% alpha and beta are derived from: input value of eps (eps = abs(nbo_enrgy)); alpha =
% sqrt(2*eps); and beta  = ((Z-Ne+1)/alpha - 1), as computed above).

B0  = intgrl_B0_num(A0, Z, alpha, Ne);

% initial parameter guess for fmincon
ini_B0    = B0;
ini_alpha = alpha;
ini_C0    = input_bcwgt(1,1);
ini_mm    = input_bcwgt(2,1);
ini_gamma = input_bcwgt(3,1);

% Nominal value for beta is beta  = ((Z-Ne+1)/alpha - 1)
if (Z == Ne) || (Z > Ne) 	% for neutral atoms or cations
	ini_beta = beta;
elseif (Z < Ne) 			% for singly-charged anions
	ini_beta = 1/alpha - 1;
end

% For documentation; alpha and beta are fit independently of eps
ini_eps   = eps;

% Set boundary conditions from input file input_bcwgt
A0LB    = input_bcwgt(4,1);
A0UB    = input_bcwgt(5,1);

B0LB 	= input_bcwgt(6,1);
B0UB 	= input_bcwgt(7,1);

C0LB 	= input_bcwgt(8,1);
C0UB 	= input_bcwgt(9,1);

mmLB 	= input_bcwgt(10,1);  
mmUB 	= input_bcwgt(11,1);  

alphaLB = input_bcwgt(12,1);
alphaUB = input_bcwgt(13,1);

betaLB  = input_bcwgt(14,1);
betaUB  = input_bcwgt(15,1);

gammaLB = input_bcwgt(16,1);
gammaUB = input_bcwgt(17,1);

% Weights for objective function
w1		= input_bcwgt(18,1);
w2		= input_bcwgt(19,1);
w3		= input_bcwgt(20,1);
w4		= input_bcwgt(21,1);
w5		= input_bcwgt(22,1);

weightSum = w1 + w2 + w3 + w4 + w5;

disp(['The sum of weights is: ',num2str(weightSum)])
errMsg2  = 'ERROR: The weights must sum to 1!';

if (abs(weightSum - 1) >= 1e-4)
    error(errMsg2)
end

% Cutoffs for fitting ranges
r1		= input_bcwgt(23,1);
r2		= input_bcwgt(24,1);

% Cutoff for long-range contribution, to avoid blow-up at nucleus
rCut    = input_bcwgt(25,1);

% --------Define objective function for fmincon as a function of x----------------

fun3 = @(x)Model3ObjFn(x, inR, inRho, Z, Ne, w1, w2, w3, w4, w5, r1, r2, rCut);

x0  = [ini_B0, ini_alpha, ini_C0, ini_mm, ini_gamma, ini_A0, ini_beta];
A   = [];
b   = [];
Aeq = [];
beq = [];

fprintf(fid1,'Small number = %f.\n\n',1.e-6);
fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha computed from input eps = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);
fprintf(fid1,'Upper and lower bound constraints for B0:    [%f,%f]\n',B0LB, B0UB);
fprintf(fid1,'Upper and lower bound constraints for alpha: [%f,%f]\n',alphaLB, alphaUB);

lb = [B0LB, alphaLB, C0LB, mmLB, gammaLB, A0LB, betaLB];
ub = [B0UB, alphaUB, C0UB, mmUB, gammaUB, A0UB, betaUB];

% Optimize parameters using fmincon
[bestx, ObjFn,exitflag,output] = fmincon(fun3,x0,A,b,Aeq,beq,lb,ub)

B0         = bestx(1);
alpha      = bestx(2);
C0         = bestx(3);
mm         = bestx(4);
gamma      = bestx(5);
A0         = bestx(6);
beta       = bestx(7);

% Exit if error in fmincon
ObjFn(isempty(ObjFn)) = "Exiting due to infeasibility";

% check for successful file close
assert(fclose(fid1)==0);

%% ------------------------------ Model 4 -------------------------------------
% Same as Model 3, with addition of Zeff as fitting parameter, in the spirit of Slater-type orbitals.
% (see GAK dissertation text for details).

elseif(model == 4)

% Import BCs and weights
input_bcwgt = readtable(input_bcwgt, 'HeaderLines',0);

% Convert to output type
input_bcwgt = table2array(input_bcwgt(:,2));

% Run analytic Model 1 to compute initial guess for B0.  Initial guesses for
% alpha and beta are derived from: input value of eps (eps = abs(nbo_enrgy)); alpha =
% sqrt(2*eps); and beta  = ((Z-Ne+1)/alpha - 1), as computed above).
B0  = intgrl_B0_num(A0, Z, alpha, Ne);

% initial parameter guess for fmincon
ini_B0    = B0;
ini_alpha = alpha;
ini_C0    = input_bcwgt(1,1);
ini_mm    = input_bcwgt(2,1);
ini_gamma = input_bcwgt(3,1);

% Nominal value for beta is beta  = ((Z-Ne+1)/alpha - 1)
if (Z == Ne) || (Z > Ne) 	  % for neutral atoms or cations
	ini_beta = beta;
elseif (Z < Ne) 			  % for anions
	ini_beta = 1/alpha - 1;
end

% For documentation; alpha and beta are fit independently of eps
ini_eps  = eps;

% Set boundary conditions from input file input_bcwgt
A0LB    = input_bcwgt(4,1);
A0UB    = input_bcwgt(5,1);

B0LB    = input_bcwgt(6,1);
B0UB    = input_bcwgt(7,1);

C0LB    = input_bcwgt(8,1);
C0UB    = input_bcwgt(9,1);

mmLB    = input_bcwgt(10,1);
mmUB    = input_bcwgt(11,1);

alphaLB = input_bcwgt(12,1);
alphaUB = input_bcwgt(13,1);

betaLB  = input_bcwgt(14,1);
betaUB  = input_bcwgt(15,1);

gammaLB = input_bcwgt(16,1);
gammaUB = input_bcwgt(17,1);

% Weights for objective function
w1      = input_bcwgt(18,1);
w2      = input_bcwgt(19,1);
w3      = input_bcwgt(20,1);
w4      = input_bcwgt(21,1);
w5      = input_bcwgt(22,1);

weightSum = w1 + w2 + w3 + w4 + w5;

disp(['The sum of weights is:',num2str(weightSum)])
errMsg2  = 'ERROR: The weights must sum to 1!';

if (abs(weightSum - 1) >= 1e-4)
    error(errMsg2)
end

% Cutoffs for fitting ranges
r1       = input_bcwgt(23,1);
r2       = input_bcwgt(24,1);

% Cutoff for long-range contribution, to avoid blow-up at nucleus
rCut     = input_bcwgt(25,1);

% New parameters for Model 4
ini_Zeff = Z;
ZeffLB   = input_bcwgt(26,1);
ZeffUB   = input_bcwgt(27,1);

% -----Define objective function for fminsearch as a function of x------

fun4 = @(x)Model4ObjFn(x, inR, inRho, Ne, w1, w2, w3, w4, w5, r1, r2, rCut);

x0 = [ini_B0, ini_alpha, ini_C0, ini_mm, ini_gamma, ini_A0, ini_beta, ini_Zeff];
A = [];
b = [];
Aeq = [];
beq = [];

fprintf(fid1,'Small number = %f.\n\n',1.e-6);
fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha from Model 1 = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);
fprintf(fid1,'Upper and lower bound constraints for B0:    [%f,%f]\n',B0LB, B0UB);
fprintf(fid1,'Upper and lower bound constraints for alpha: [%f,%f]\n',alphaLB, alphaUB);

lb = [B0LB, alphaLB, C0LB, mmLB, gammaLB, A0LB, betaLB, ZeffLB];
ub = [B0UB, alphaUB, C0UB, mmUB, gammaUB, A0UB, betaUB, ZeffUB];

% Optimize B0 and alpha parameters using fminsearch
% options  = optimset('Display','iter','PlotFcns',@optimplotfval,'TolX',1e-6);
% bestx    = fminsearch(fun,x0,options);
[bestx, ObjFn,exitflag,output] = fmincon(fun4,x0,A,b,Aeq,beq,lb,ub)

B0         = bestx(1);
alpha      = bestx(2);
C0         = bestx(3);
mm         = bestx(4);
gamma      = bestx(5);
A0         = bestx(6);
beta       = bestx(7);
Zeff       = bestx(8);

betaFFAZ    = ((Zeff-Ne+1)/alpha - 1);    % beta computed From Fitted Alpha and Zeff

% exit if encounter empty objextive function
ObjFn(isempty(ObjFn)) = "Exiting due to infeasibility";

fprintf(fid1,'Optimized B0 from Model 4                 = %f\n',B0);
fprintf(fid1,'Optimized Zeff from Model 4               = %f\n',Zeff);
fprintf(fid1,'Optimized alpha from Model 4              = %f\n',alpha);
fprintf(fid1,'Optimized beta from Model 4               = %f\n',beta);
fprintf(fid1,'Beta computed from Model 4 alpha and Zeff = %f\n',betaFFAZ);
fprintf(fid1,'Optimized C0 from Model 4                 = %f\n',C0);
fprintf(fid1,'Optimized mm from Model 4                 = %f\n',mm);
fprintf(fid1,'Optimized gamma from Model 4              = %f\n',gamma);
assert(fclose(fid1)==0);   % check for successful file close

%% ------------------------------ Model 5 -------------------------------------
% Same as Model 4, with addition of a 2p radial orbital density, and
% new fitting parameters D0 and eta (see GAK dissertation text for details).

elseif(model == 5)

% Import BCs and weights
input_bcwgt = readtable(input_bcwgt, 'HeaderLines',0);

% Convert to output type
input_bcwgt = table2array(input_bcwgt(:,2))

% Run analytic Model 1 to compute initial guess for B0 and alpha parameter.
B0  = intgrl_B0_num(A0, Z, alpha, Ne);

% initial parameter guess for fmincon
% Note: ini_A0 is set above, when g09 or dbar data is read in. 

ini_B0    = B0;
ini_alpha = alpha;
ini_C0    = input_bcwgt(1,1);
ini_mm    = input_bcwgt(2,1);
ini_gamma = input_bcwgt(3,1);

% For documentation; alpha and beta are fit independently of eps
ini_eps   = eps;                

% Nominal value for beta is beta  = ((Z-Ne+1)/alpha - 1)
% if (Z == Ne) || (Z > Ne) 	  % for neutral atoms or cations
% 	ini_beta = beta;
% elseif (Z < Ne) 			  % for anions
% 	ini_beta = 1/alpha - 1;
% end

ini_beta = beta;

% Set boundary conditions from input file input_bcwgt
A0LB    = input_bcwgt(4,1);
A0UB    = input_bcwgt(5,1);

B0LB    = input_bcwgt(6,1);
B0UB    = input_bcwgt(7,1);

C0LB    = input_bcwgt(8,1);
C0UB    = input_bcwgt(9,1);

mmLB    = input_bcwgt(10,1);
mmUB    = input_bcwgt(11,1);  

alphaLB = input_bcwgt(12,1);
alphaUB = input_bcwgt(13,1);

betaLB  = input_bcwgt(14,1);
betaUB  = input_bcwgt(15,1);

gammaLB = input_bcwgt(16,1);
gammaUB = input_bcwgt(17,1);

% Weights for objective function
w1      = input_bcwgt(18,1);
w2      = input_bcwgt(19,1);
w3      = input_bcwgt(20,1);
w4      = input_bcwgt(21,1);
w5      = input_bcwgt(22,1);

weightSum = w1 + w2 + w3 + w4 + w5;

disp(['The sum of weights is:',num2str(weightSum)])
errMsg2  = 'ERROR: The weights must sum to 1!';

if (abs(weightSum - 1) >= 1e-4)
    error(errMsg2)
end

% Cutoffs for ranges
r1        = input_bcwgt(23,1);
r2        = input_bcwgt(24,1);
rCut      = input_bcwgt(25,1);

% Zeff
ini_Zeff  = Z;
ZeffLB    = input_bcwgt(26,1);
ZeffUB    = input_bcwgt(27,1);

% New parameters for Model 5
ini_D0    = input_bcwgt(28,1);
ini_eta   = input_bcwgt(29,1);

etaLB     = input_bcwgt(30,1);
etaUB     = input_bcwgt(31,1);

D0LB      = input_bcwgt(32,1);
D0UB      = input_bcwgt(33,1);

% --------Define objective function for fmincon as a function of x----------------

fun5 = @(x)Model5ObjFn(x, inR, inRho, Ne, w1, w2, w3, w4, w5, r1, r2, rCut);

x0 = [ini_B0, ini_alpha, ini_C0, ini_mm, ini_gamma, ini_A0, ini_beta, ini_Zeff, ini_D0, ini_eta];
A = [];
b = [];
Aeq = [];
beq = [];

fprintf(fid1,'Small number = %f.\n\n',1.e-6);
fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha from Model 1 = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);

lb = [B0LB, alphaLB, C0LB, mmLB, gammaLB, A0LB, betaLB, ZeffLB, D0LB, etaLB]
ub = [B0UB, alphaUB, C0UB, mmUB, gammaUB, A0UB, betaUB, ZeffUB, D0UB, etaUB]

% Constrained optimization of model parameters using Matlab function fmincon
%[bestx, ObjFn,exitflag,output] = fmincon(fun5,x0,A,b,Aeq,beq,lb,ub)

%opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
%[x_global, ~] = particleswarm(fun5, length(x0), lb, ub, opts1);

%opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
%[bestx, ObjFn,exitflag,output] = fmincon(fun5, x_global, A, b, Aeq, beq, lb, ub, [], opts2);

if useGlobalSearch
    % Step 1: Global search
    opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
    [x_global, ~] = particleswarm(fun5, length(x0), lb, ub, opts1);

    % === Log particle swarm starting point ===
    startParamFile = sprintf('%sstartParams_model%s_%s.txt', plot_path, modelNam, dateStrTime);
    fidSP = fopen(startParamFile, 'wt');
    paramNames = ["B0","alpha","C0","mm","gamma","A0","beta","Zeff","D0","eta"];
    for k = 1:length(x_global)
        fprintf(fidSP, '%s\t%.8g\n', paramNames(k), x_global(k));
    end
    fprintf(fidSP, '\n--- Bounds and Weights ---\n');
    fullNames  = ["B0LB","B0UB","alphaLB","alphaUB","C0LB","C0UB","mmLB","mmUB","betaLB","betaUB","gammaLB","gammaUB",...
                  "A0LB","A0UB","D0LB","D0UB","etaLB","etaUB","ZeffLB","ZeffUB","w1","w2","w3","w4","w5","r1","r2","rCut"];
    fullValues = [B0LB,B0UB,alphaLB,alphaUB,C0LB,C0UB,mmLB,mmUB,betaLB,betaUB,gammaLB,gammaUB,...
                  A0LB,A0UB,D0LB,D0UB,etaLB,etaUB,ZeffLB,ZeffUB,w1,w2,w3,w4,w5,r1,r2,rCut];
    for k = 1:length(fullNames)
        fprintf(fidSP, '%s\t%.8g\n', fullNames(k), fullValues(k));
    end
    fclose(fidSP);

    % Step 2: Local refinement
    opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',1e5, ...    % allow 100,000 evals
    'MaxIterations',5000);               % allow more iterations
    [bestx, ObjFn,exitflag,output] = fmincon(fun5, x_global, A, b, Aeq, beq, lb, ub, [], opts2);
else
    [bestx, ObjFn,exitflag,output] = fmincon(fun5, x0, A, b, Aeq, beq, lb, ub);
end

B0         = bestx(1);
alpha      = bestx(2);
C0         = bestx(3);
mm         = bestx(4);
gamma      = bestx(5);
A0         = bestx(6);
beta       = bestx(7);
Zeff       = bestx(8);
D0         = bestx(9);
eta        = bestx(10);

% Variable name is an acronym: beta computed From Fitted Alpha and Zeff
betaFFAZ   = ((Zeff-Ne+1)/alpha - 1);

% exit if encounter empty objective function
ObjFn(isempty(ObjFn)) = "Exiting due to infeasibility";

fprintf(fid1,'Optimized A0 from Model 5                 = %f\n',A0);
fprintf(fid1,'Optimized B0 from Model 5                 = %f\n',B0);
fprintf(fid1,'Optimized C0 from Model 5                 = %f\n',C0);
fprintf(fid1,'Optimized D0 from Model 5                 = %f\n',D0);
fprintf(fid1,'Optimized Zeff from Model 5               = %f\n',Zeff);
fprintf(fid1,'Optimized alpha from Model 5              = %f\n',alpha);
fprintf(fid1,'Optimized beta from Model 5               = %f\n',beta);
fprintf(fid1,'Beta computed from Model 5 alpha and Zeff = %f\n',betaFFAZ);
fprintf(fid1,'Optimized mm from Model 5                 = %f\n',mm);
fprintf(fid1,'Optimized gamma from Model 5              = %f\n',gamma);
assert(fclose(fid1)==0);   % check for successful file close

%% ------------------------------ Model 6 -------------------------------------
% Same as Model 5, with addition of linear fit for beta, and new parameter bb
% (see GAK dissertation text for details).

elseif(model == 6)

% Import BCs and weights
input_bcwgt = readtable(input_bcwgt, 'HeaderLines',0);

% Convert to output type
input_bcwgt = table2array(input_bcwgt(:,2))

% Run analytic Model 1 to compute initial guess for B0 and alpha parameter.
B0  = intgrl_B0_num(A0, Z, alpha, Ne);

% initial parameter guess for fmincon
% Note: ini_A0 is set above, when g09 or dbar data is read in. 

ini_B0    = B0;
ini_alpha = alpha;
ini_C0    = input_bcwgt(1,1);
ini_mm    = input_bcwgt(2,1);
ini_gamma = input_bcwgt(3,1);

% For documentation only; alpha and beta are fit independently of eps
ini_eps   = eps;                

% Set boundary conditions from input file input_bcwgt
A0LB    = input_bcwgt(4,1);
A0UB    = input_bcwgt(5,1);

B0LB    = input_bcwgt(6,1);
B0UB    = input_bcwgt(7,1);

C0LB    = input_bcwgt(8,1);
C0UB    = input_bcwgt(9,1);

mmLB    = input_bcwgt(10,1);
mmUB    = input_bcwgt(11,1);  

alphaLB = input_bcwgt(12,1);
alphaUB = input_bcwgt(13,1);

betaLB  = input_bcwgt(14,1);
betaUB  = input_bcwgt(15,1);

gammaLB = input_bcwgt(16,1);
gammaUB = input_bcwgt(17,1);

% Nominal value for beta is beta  = ((Z-Ne+1)/alpha - 1)
%
% Initial attempt:
%    if (Z == Ne) || (Z > Ne) 	  % for neutral atoms or cations
% 	   ini_beta = beta;
%    elseif (Z < Ne) 			  % for anions; try initial value .ne. nominal value of -1
% 	   ini_beta = 1/alpha - 1;
%    end
%
% Instead, the following enables fixing beta = betaLB = betaUB = value for corresponding *neutral* atom.

ini_beta = betaLB;

% Weights for objective function
w1      = input_bcwgt(18,1);
w2      = input_bcwgt(19,1);
w3      = input_bcwgt(20,1);
w4      = input_bcwgt(21,1);
w5      = input_bcwgt(22,1);

weightSum = w1 + w2 + w3 + w4 + w5;

disp(['The sum of weights is:',num2str(weightSum)])
errMsg2  = 'ERROR: The weights must sum to 1!';

if (abs(weightSum - 1) >= 1e-4)
    error(errMsg2)
end

% Cutoffs for ranges
r1        = input_bcwgt(23,1);
r2        = input_bcwgt(24,1);
rCut      = input_bcwgt(25,1);

% Zeff
ini_Zeff  = Z;
ZeffLB    = input_bcwgt(26,1);
ZeffUB    = input_bcwgt(27,1);

% New parameters for Model 5
ini_D0    = input_bcwgt(28,1);
ini_eta   = input_bcwgt(29,1);

etaLB     = input_bcwgt(30,1);
etaUB     = input_bcwgt(31,1);

D0LB      = input_bcwgt(32,1);
D0UB      = input_bcwgt(33,1);

% New parameter for Model 6
ini_bb    = input_bcwgt(34,1);

bbLB      = input_bcwgt(35,1);
bbUB      = input_bcwgt(36,1);

% --------Define objective function for fmincon as a function of x----------------

fun6 = @(x)Model6ObjFn(x, inR, inRho, Ne, w1, w2, w3, w4, w5, r1, r2, rCut);

x0 = [ini_B0, ini_alpha, ini_C0, ini_mm, ini_gamma, ini_A0, ini_beta, ini_Zeff, ini_D0, ini_eta, ini_bb];
A = [];
b = [];
Aeq = [];
beq = [];

fprintf(fid1,'Small number = %f.\n\n',1.e-6);
fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha from Model 1 = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);

lb = [B0LB, alphaLB, C0LB, mmLB, gammaLB, A0LB, betaLB, ZeffLB, D0LB, etaLB, bbLB]
ub = [B0UB, alphaUB, C0UB, mmUB, gammaUB, A0UB, betaUB, ZeffUB, D0UB, etaUB, bbUB]

% Constrained optimization of model parameters using Matlab function fmincon
% [bestx, ObjFn,exitflag,output] = fmincon(fun6,x0,A,b,Aeq,beq,lb,ub)

% Step 1: Global search
%opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
%[x_global, ~] = particleswarm(fun6, length(x0), lb, ub, opts1);

% Step 2: Refine with fmincon
%opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
%[bestx, ObjFn,exitflag,output] = fmincon(fun6, x_global, A, b, Aeq, beq, lb, ub, [], opts2);


if useGlobalSearch
    opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
    [x_global, ~] = particleswarm(fun6, length(x0), lb, ub, opts1);

    % === Log particle swarm starting point ===
    startParamFile = sprintf('%sstartParams_model%s_%s.txt', plot_path, modelNam, dateStrTime);
    fidSP=fopen(startParamFile,'wt');
    paramNames=["B0","alpha","C0","mm","gamma","A0","beta","Zeff","D0","eta","bb"];
    for k=1:length(x_global)
        fprintf(fidSP,'%s\t%.8g\n',paramNames(k),x_global(k));
    end
    fprintf(fidSP,'\n--- Bounds and Weights ---\n');
    fullNames=["B0LB","B0UB","alphaLB","alphaUB","C0LB","C0UB","mmLB","mmUB","betaLB","betaUB","gammaLB","gammaUB",...
               "A0LB","A0UB","D0LB","D0UB","etaLB","etaUB","bbLB","bbUB","ZeffLB","ZeffUB","w1","w2","w3","w4","w5","r1","r2","rCut"];
    fullValues=[B0LB,B0UB,alphaLB,alphaUB,C0LB,C0UB,mmLB,mmUB,betaLB,betaUB,gammaLB,gammaUB,...
                A0LB,A0UB,D0LB,D0UB,etaLB,etaUB,bbLB,bbUB,ZeffLB,ZeffUB,w1,w2,w3,w4,w5,r1,r2,rCut];
    for k=1:length(fullNames)
        fprintf(fidSP,'%s\t%.8g\n',fullNames(k),fullValues(k));
    end
    fclose(fidSP);

    opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
    [bestx, ObjFn,exitflag,output] = fmincon(fun6, x_global, A, b, Aeq, beq, lb, ub, [], opts2);
else
    [bestx, ObjFn,exitflag,output] = fmincon(fun6,x0,A,b,Aeq,beq,lb,ub)
end

B0         = bestx(1);
alpha      = bestx(2);
C0         = bestx(3);
mm         = bestx(4);
gamma      = bestx(5);
A0         = bestx(6);
beta       = bestx(7);
Zeff       = bestx(8);
D0         = bestx(9);
eta        = bestx(10);
bb         = bestx(11);

% Variable name is an acronym: beta computed From Fitted Alpha and Zeff
betaFFAZ   = ((Zeff-Ne+1)/alpha - 1);

% exit if encounter empty objective function
ObjFn(isempty(ObjFn)) = "Exiting due to infeasibility";

fprintf(fid1,'Optimized B0 from Model 6                 = %f\n',B0);
fprintf(fid1,'Optimized Zeff from Model 6               = %f\n',Zeff);
fprintf(fid1,'Optimized alpha from Model 6              = %f\n',alpha);
fprintf(fid1,'Optimized beta from Model 6               = %f\n',beta);
fprintf(fid1,'Beta computed from Model 6 alpha and Zeff = %f\n',betaFFAZ);
fprintf(fid1,'Optimized C0 from Model 6                 = %f\n',C0);
fprintf(fid1,'Optimized mm from Model 6                 = %f\n',mm);
fprintf(fid1,'Optimized gamma from Model 6              = %f\n',gamma);
fprintf(fid1,'Optimized D0 from Model 6                 = %f\n',D0);
fprintf(fid1,'Optimized eta from Model 6                = %f\n',eta);
fprintf(fid1,'Optimized bb from Model 6                 = %f\n',bb);

% Print short-range and ending beta points from the fit, for tabulation and
% to illustrate transition from neutral atom effective beta to purely anionic beta (-1).

betav   = (inR==0)*1 + (inR>0)*beta;
betaLin = betav + bb.*inR;
lbl     = length(betaLin);

fprintf(fid1,'\n');
fprintf(fid1,'beta(0)= %f.\n',beta);
fprintf(fid1,'beta(%f)= %f.\n',inR(lbl),betaLin(lbl));

assert(fclose(fid1)==0);   % check for successful file close
%% ------------------------------ Model 7 -------------------------------------
% Same as Model 6, with addition of 3s state, and new parameters E0 and
% delta, pp, qq, tt
% (see GAK dissertation text for details).

elseif(model == 7)

% Import BCs and weights
input_bcwgt = readtable(input_bcwgt, 'HeaderLines',0);

% Convert to output type
input_bcwgt = table2array(input_bcwgt(:,2))

% Run analytic Model 1 to compute initial guess for B0 and alpha parameter.
B0  = intgrl_B0_num(A0, Z, alpha, Ne);

% initial parameter guess for fmincon
% Note: ini_A0 is set above, when g09 or dbar data is read in. 

ini_B0    = B0;
ini_alpha = alpha;
ini_C0    = input_bcwgt(1,1);
ini_mm    = input_bcwgt(2,1);
ini_gamma = input_bcwgt(3,1);

% For documentation only; alpha and beta are fit independently of eps
ini_eps   = eps;                

% Set boundary conditions from input file input_bcwgt
A0LB    = input_bcwgt(4,1);
A0UB    = input_bcwgt(5,1);

B0LB    = input_bcwgt(6,1);
B0UB    = input_bcwgt(7,1);

C0LB    = input_bcwgt(8,1);
C0UB    = input_bcwgt(9,1);

mmLB    = input_bcwgt(10,1);
mmUB    = input_bcwgt(11,1);  

alphaLB = input_bcwgt(12,1);
alphaUB = input_bcwgt(13,1);

betaLB  = input_bcwgt(14,1);
betaUB  = input_bcwgt(15,1);

gammaLB = input_bcwgt(16,1);
gammaUB = input_bcwgt(17,1);

% Nominal value for beta is beta  = ((Z-Ne+1)/alpha - 1)
%
% Initial attempt:
%    if (Z == Ne) || (Z > Ne) 	  % for neutral atoms or cations
% 	   ini_beta = beta;
%    elseif (Z < Ne) 			  % for anions; try initial value .ne. nominal value of -1
% 	   ini_beta = 1/alpha - 1;
%    end
%
% Instead, the following enables fixing beta = betaLB = betaUB = value for corresponding *neutral* atom.

ini_beta = betaLB;

% Weights for objective function
w1      = input_bcwgt(18,1);
w2      = input_bcwgt(19,1);
w3      = input_bcwgt(20,1);
w4      = input_bcwgt(21,1);
w5      = input_bcwgt(22,1);

weightSum = w1 + w2 + w3 + w4 + w5;

disp(['The sum of weights is:',num2str(weightSum)])
errMsg2  = 'ERROR: The weights must sum to 1!';

if (abs(weightSum - 1) >= 1e-4)
    error(errMsg2)
end

% Cutoffs for ranges
r1        = input_bcwgt(23,1);
r2        = input_bcwgt(24,1);
rCut      = input_bcwgt(25,1);

% Zeff
ini_Zeff  = Z;
ZeffLB    = input_bcwgt(26,1);
ZeffUB    = input_bcwgt(27,1);

% New parameters for Model 5
ini_D0    = input_bcwgt(28,1);
ini_eta   = input_bcwgt(29,1);

etaLB     = input_bcwgt(30,1);
etaUB     = input_bcwgt(31,1);

D0LB      = input_bcwgt(32,1);
D0UB      = input_bcwgt(33,1);

ini_bb    = input_bcwgt(34,1);

bbLB      = input_bcwgt(35,1);
bbUB      = input_bcwgt(36,1);

% New parameters for Model 7
ini_E0    = input_bcwgt(37,1);
E0LB      = input_bcwgt(38,1);
E0UB      = input_bcwgt(39,1);

ini_delta    = input_bcwgt(40,1);
deltaLB      = input_bcwgt(41,1);
deltaUB      = input_bcwgt(42,1);

ini_pp    = input_bcwgt(43,1);
ppLB      = input_bcwgt(44,1);
ppUB      = input_bcwgt(45,1);

ini_qq    = input_bcwgt(46,1);
qqLB      = input_bcwgt(47,1);
qqUB      = input_bcwgt(48,1);

ini_tt    = input_bcwgt(49,1);
ttLB      = input_bcwgt(50,1);
ttUB      = input_bcwgt(51,1);



% --------Define objective function for fmincon as a function of x----------------

fun7 = @(x)Model7ObjFn(x, inR, inRho, Ne, w1, w2, w3, w4, w5, r1, r2, rCut);

x0 = [ini_B0, ini_alpha, ini_C0, ini_mm, ini_gamma, ini_A0, ini_beta, ini_Zeff, ini_D0, ini_eta, ini_bb, ini_E0, ini_delta, ini_pp, ini_qq, ini_tt];
A = [];
b = [];
Aeq = [];
beq = [];

fprintf(fid1,'Small number = %f.\n\n',1.e-6);
fprintf(fid1,'Initial guess for B0 from Model 1 = %f.\n',B0);
fprintf(fid1,'Initial guess for alpha from Model 1 = %f.\n',alpha);
fprintf(fid1,'Initial guess for beta computed from Model 1 alpha = %f.\n\n',beta);

lb = [B0LB, alphaLB, C0LB, mmLB, gammaLB, A0LB, betaLB, ZeffLB, D0LB, etaLB, bbLB, E0LB, deltaLB, ppLB, qqLB, ttLB]
ub = [B0UB, alphaUB, C0UB, mmUB, gammaUB, A0UB, betaUB, ZeffUB, D0UB, etaUB, bbUB, E0UB, deltaUB, ppUB, qqUB, ttUB]

% Constrained optimization of model parameters using Matlab function fmincon
%[bestx, ObjFn,exitflag,output] = fmincon(fun7,x0,A,b,Aeq,beq,lb,ub)

%inclusion of initial global search, then local search to handle larger
%parameter space
% Step 1: Global search
%opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
%[x_global, ~] = particleswarm(fun7, length(x0), lb, ub, opts1);

% Step 2: Refine with fmincon
%opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
%[bestx, ObjFn,exitflag,output] = fmincon(fun7, x_global, A, b, Aeq, beq, lb, ub, [], opts2);

if useGlobalSearch
    % Step 1: Global search
    opts1 = optimoptions('particleswarm','Display','iter','SwarmSize',200);
    [x_global, ~] = particleswarm(fun7, length(x0), lb, ub, opts1);

    % === Log particle swarm starting point ===
    startParamFile = sprintf('%sstartParams_model%s_%s.txt', plot_path, modelNam, dateStrTime);
    fidSP = fopen(startParamFile,'wt');

    paramNames = ["B0","alpha","C0","mm","gamma","A0","beta","Zeff","D0","eta","bb","E0","delta","pp","qq","tt"];
    for k = 1:length(x_global)
        fprintf(fidSP,'%s\t%.8g\n',paramNames(k),x_global(k));
    end

    fprintf(fidSP,'\n--- Bounds and Weights ---\n');
    fullNames = ["C0","mm","gamma","A0LB","A0UB","B0LB","B0UB","C0LB","C0UB","mmLB","mmUB","alphaLB","alphaUB","betaLB","betaUB","gammaLB","gammaUB",...
                 "w1","w2","w3","w4","w5","r1","r2","rCut","ZeffLB","ZeffUB","D0","eta","etaLB","etaUB","D0LB","D0UB",...
                 "bb","bbLB","bbUB","E0","E0LB","E0UB","delta","deltaLB","deltaUB","pp","ppLB","ppUB","qq","qqLB","qqUB","tt","ttLB","ttUB"];
    fullValues = [ini_C0,ini_mm,ini_gamma,A0LB,A0UB,B0LB,B0UB,C0LB,C0UB,mmLB,mmUB,alphaLB,alphaUB,betaLB,betaUB,gammaLB,gammaUB, ...
                  w1,w2,w3,w4,w5,r1,r2,rCut,ZeffLB,ZeffUB,ini_D0,ini_eta,etaLB,etaUB,D0LB,D0UB, ...
                  ini_bb,bbLB,bbUB,ini_E0,E0LB,E0UB,ini_delta,deltaLB,deltaUB,ini_pp,ppLB,ppUB,ini_qq,qqLB,qqUB,ini_tt,ttLB,ttUB];

    for k = 1:length(fullNames)
        fprintf(fidSP,'%s\t%.8g\n',fullNames(k),fullValues(k));
    end

    fclose(fidSP);

    % Step 2: Local refinement
    opts2 = optimoptions('fmincon','Display','iter','Algorithm','interior-point', ...
    'MaxFunctionEvaluations',1e5, ...    % allow 100,000 evals
    'MaxIterations',5000);               % allow more iterations
    [bestx, ObjFn,exitflag,output] = fmincon(fun7, x_global, A, b, Aeq, beq, lb, ub, [], opts2);
else
    % Skip global search
    [bestx, ObjFn,exitflag,output] = fmincon(fun7,x0,A,b,Aeq,beq,lb,ub);
end

B0         = bestx(1);
alpha      = bestx(2);
C0         = bestx(3);
mm         = bestx(4);
gamma      = bestx(5);
A0         = bestx(6);
beta       = bestx(7);
Zeff       = bestx(8);
D0         = bestx(9);
eta        = bestx(10);
bb         = bestx(11);
E0         = bestx(12);
delta      = bestx(13);
pp         = bestx(14);
qq         = bestx(15);
tt         = bestx(16);

% Variable name is an acronym: beta computed From Fitted Alpha and Zeff
betaFFAZ   = ((Zeff-Ne+1)/alpha - 1);

% exit if encounter empty objective function
ObjFn(isempty(ObjFn)) = "Exiting due to infeasibility";

fprintf(fid1,'Optimized B0 from Model 7                 = %f\n',B0);
fprintf(fid1,'Optimized Zeff from Model 7               = %f\n',Zeff);
fprintf(fid1,'Optimized alpha from Model 7              = %f\n',alpha);
fprintf(fid1,'Optimized beta from Model 7               = %f\n',beta);
fprintf(fid1,'Beta computed from Model 7 alpha and Zeff = %f\n',betaFFAZ);
fprintf(fid1,'Optimized C0 from Model 7                 = %f\n',C0);
fprintf(fid1,'Optimized mm from Model 7                 = %f\n',mm);
fprintf(fid1,'Optimized gamma from Model 7              = %f\n',gamma);
fprintf(fid1,'Optimized D0 from Model 7                 = %f\n',D0);
fprintf(fid1,'Optimized eta from Model 7                = %f\n',eta);
fprintf(fid1,'Optimized bb from Model 7                 = %f\n',bb);
fprintf(fid1,'Optimized E0 from Model 7                 = %f\n',E0);
fprintf(fid1,'Optimized delta from Model 7              = %f\n',delta);
fprintf(fid1,'Optimized pp from Model 7                 = %f\n',pp);
fprintf(fid1,'Optimized qq from Model 7                 = %f\n',qq);
fprintf(fid1,'Optimized tt from Model 7                 = %f\n',tt);

% Print short-range and ending beta points from the fit, for tabulation and
% to illustrate transition from neutral atom effective beta to purely anionic beta (-1).

betav   = (inR==0)*1 + (inR>0)*beta;
betaLin = betav + bb.*inR;
lbl     = length(betaLin);

fprintf(fid1,'\n');
fprintf(fid1,'beta(0)= %f.\n',beta);
fprintf(fid1,'beta(%f)= %f.\n',inR(lbl),betaLin(lbl));

assert(fclose(fid1)==0);   % check for successful file close
end
% -------------------------------- End model fitting --------------------------

% ---------------------------- Begin model evaluations ------------------------

%% Compute the model density from returned parameters to check the quality of the fit, 
%  To check the quality of the fit, compute the model density from the  returned parameters

if(model == 1) || (model == 2)
    model_fit  = density_model_1_2_eval_fn(A0, B0, Z, alpha, Ne, inR);

    radial_dBar      = rdBar.^2.*rhodBar;                   % RDF for dBar
    radial_model_fit = fourPi*inR.^2.*model_fit;            % RDF for analytic model with fitted B0
    
    % Calculate the integrated number of electrons for all three data files
    % Note that except for spherically-symmetric atoms, Int_g09 will not integrate perfectly to total #e-
    
    Int_dBar   = double(trapz(rdBar, radial_dBar));          % #e dBar
    Int_model  = double(trapz(inR, radial_model_fit));       % #e analytic model fit to dBar
    
elseif(model >= 3) 
  if(model == 4) || (model == 5) || (model == 6) || (model == 7)
    Z = Zeff; 
  end
    
  if(model == 3) || (model == 4)
     model_fit = density_model_3_4_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, inR, rCut);
  elseif(model == 5)
     model_fit = density_model_5_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, D0, eta, inR, rCut);
  elseif(model == 6)
     model_fit = density_model_6_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, D0, eta, bb, inR, rCut);
  elseif(model == 7)
     model_fit = density_model_7_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, D0, eta, bb, E0, delta, pp, qq, tt, inR, rCut);
  end

% Generate radial distribution functions 4*pi*r^2*\rho(r)
% Note that dBar as defined already includes factor of 4*pi

radial_dBar      = rdBar.^2.*rhodBar;                   % RDF for dBar
radial_model_fit = fourPi*inR.^2.*model_fit;            % RDF for analytic model with fitted B0

% Calculate the integrated number of electrons for all three data files
% Note that except for spherically-symmetric atoms, Int_g09 will not integrate perfectly to total #e-

Int_dBar   = double(trapz(rdBar, radial_dBar));          % #e dBar
Int_model  = double(trapz(inR, radial_model_fit));       % #e analytic model fit to dBar

end
% ----------------------------- End model evaluations -------------------------

% Generate plots of model density, ln density, and RDF; compare with input
% density (G'09 zcut or dBar) 
close all;
figure('Name','Densities');

% dBar density
rhodBarOver4pi = rhodBar./fourPi;
plot(rdBar,(rhodBarOver4pi),'-g','LineWidth',2,'MarkerIndices',1:50:length(rhodBar));
hold on

% model density
plot(inR,model_fit,'-b','LineWidth',2,'MarkerIndices',1:40:length(model_fit));
hold on

legend('dBar', ...
        'Analytical model', ...
        'Location','NorthEast')
    
% grid on
% grid minor
hold off

xlabel('Radial Distance, $\rm{r~(au)} $ ','interpreter','LaTex')
ylabel('$\ln \rho$ (au)','interpreter','LaTex')
titlename = char(sprintf('Density, %s', atomicstate));
title(titlename);
% xlim([0 3])
% ylim([0 500])

% Save as pdf file
filename = char(sprintf('density_%s_%s', modelNam, dateStrTime));

% saveas(gcf, filename1b, 'pdf'); %toggle if need to save in current directory
savePlot(plot_path,filename)

% ----------------- Plot ln of G09 zcut,dBar,and model densities ------------

figure('Name','ln densities');


% dBar density
lnrhodBarOver4pi = log(rhodBar./fourPi);
plot(rdBar,lnrhodBarOver4pi,'-g','LineWidth',2,'MarkerIndices',1:50:length(lnrhodBarOver4pi));
hold on

% model density
ln_model_fit= log(model_fit);
plot(inR,ln_model_fit,'-b','LineWidth',2,'MarkerIndices',1:40:length(ln_model_fit));
hold on

legend('dBar', ...
        'Analytical model', ...
        'Location','NorthEast')
% grid on
% grid minor
hold off

xlabel('Radial Distance, $\rm{r~(au)} $ ','interpreter','LaTex')
ylabel('$\ln \rho$ (au)','interpreter','LaTex')
titlename = char(sprintf('ln density, %s', atomicstate));
title(titlename);
% xlim([0 8])
% ylim([0 500])

% Save as pdf file
filename = char(sprintf('ln_density_%s_%s', modelNam, dateStrTime));
% saveas( gcf(), filename1b, 'pdf' ); %toggle if need to save in current directory
savePlot(plot_path,filename)

% %**************Plot and save the radial distrubution funtions.**********
rdBar        = [0; rdBar];
radial_dBar  = [0; radial_dBar];
inR         = [0; inR];
radial_model_fit = [0; radial_model_fit];

figure('Name','RDF'),

% dBar data
plot(rdBar,radial_dBar,'-g','LineWidth',2,'MarkerIndices',1:60:length(radial_dBar));
hold on

% model data
plot(inR,radial_model_fit,'-b','LineWidth',2,'MarkerIndices',1:60:length(radial_model_fit));
hold on

legend('dBar',...
       'Analytical model',...
       'Location','NorthEast');
% grid on
% grid minor
hold off

xlabel('Radial Distance, $\rm{r~(au)} $ ','interpreter','LaTex')
ylabel('$4\pi r^2 \rho(r)$ ','interpreter','LaTex')
titlename = char(sprintf('Radial Distribution Function, %s', atomicstate));
title(titlename);
% xlim([0 5])

% Save RDF plot as .pdf
filename = char(sprintf('RDF_%s_%s', modelNam, dateStrTime));
% saveas( gcf(), filename2b, 'pdf' ); %toggle if need to save in current directory
savePlot(plot_path,filename)

%% For Models 3, 4, 5, and 6 compute and plot final model component densities

if(model >= 3)
       
    if(model < 5)        % 2p contribution only in Models 5 and 6
        D0       = 0.0;                              
        eta      = 0.0;
        D0LB     = 0.0;
        D0UB     = 0.0;
        etaLB    = 0.0;
        etaUB    = 0.0;
        ini_D0   = 0.0;
        ini_eta  = 0.0;
        etaLB    = 0.0;
        etaUB    = 0.0;
        D0LB     = 0.0;
        D0UB     = 0.0;
        bestx(9) = 0.0;
        bestx(10)= 0.0;
    end
    
    if(model ~= 6 && model ~= 7)
        bb     = 0.0;
        ini_bb = 0.0;
        bbLB   = 0.0;
        bbUB   = 0.0;
        bestx(11) = 0.0;
    end

    if(model ~= 7)
        E0     = 0.0;
        ini_E0 = 0.0;
        E0LB   = 0.0;
        E0UB   = 0.0;
        bestx(12) = 0.0;
        delta     = 0.0;
        ini_delta = 0.0;
        deltaLB   = 0.0;
        deltaUB   = 0.0;
        bestx(13) = 0.0;
        pp     = 0.0;
        ini_pp = 0.0;
        ppLB   = 0.0;
        ppUB   = 0.0;
        bestx(14) = 0.0;
        qq     = 0.0;
        ini_qq = 0.0;
        qqLB   = 0.0;
        qqUB   = 0.0;
        bestx(15) = 0.0;
        tt     = 0.0;
        ini_tt = 0.0;
        ttLB   = 0.0;
        ttUB   = 0.0;
        bestx(16) = 0.0;
    end
    
    if(model == 4) || (model == 5) || (model == 6) || (model == 7)
        Z = Zeff;
    end
    
    betaFFAZ   = ((Z-Ne+1)/alpha - 1);           % beta computed From Fitted Alpha and Zeff
    
    F1         = A0.*(exp(-2*Z*inR));            % short-range cusp
    density1R  = fourPi*inR.^2.*F1;
    
    F2_a = (mm - Z*inR);
    F2_b = F2_a .* F2_a;
    F2   = C0*F2_b.*exp(-2*gamma*inR);           % 2s shell structure  
    density2R  = fourPi*inR.^2.*F2;

    F3_a   = Z*inR;
    F3_b   = F3_a .* F3_a;
    F3     = D0*F3_b.*exp(-2*eta*inR);           % 2p shell structure  
    density3R  = fourPi*inR.^2.*F3;
    

    F5_a   = pp-(qq*Z*inR);
    F5_b   = (tt*Z*inR);
    F5_c   = (F5_b).*(F5_b);
    F5_d   = F5_a + F5_c;
    F5_e   = (F5_d).*(F5_d);
    F5     = E0*F5_e.*exp(-2*delta*inR);         % 3s shell structure  
    density5R  = fourPi*inR.^2.*F5;
   
    betav   = (inR==0)*1 + (inR>0)*beta;
    if (model == 6) || (model == 7)
       betaLin = betav + bb.*inR;
       betav = betaLin;
    end 
    F4      = (inR>rCut).*(B0.*((inR).^(2.*betav)).*(exp(-2.*alpha.*(inR))));  
    density4R = fourPi*inR.^2.*F4;               % long-range decay

    totalDensityR = density1R + density2R + density3R + density4R + density5R;
    totalDensityln = log(F1+F2+F3+F4+F5);
    
% Effective beta plot
    if (model == 6) || (model == 7)
      figure('Name','betaForAnion');
      lbv    = length(betav);
      inRP   = inR;
      betavP = [beta; betav(2:lbv)];
      plot(inRP,betavP,'-b','LineWidth',2,'MarkerIndices',1:60:lbv);
      title('Effective $\beta$ for anion','interpreter','LaTex');
      xlabel('Radial Distance, $\rm{r~(au)} $ ','interpreter','LaTex')
      ylabel('$\beta(r)$ ','interpreter','LaTex');
      filename = char(sprintf('betaForAnion_%s_%s', modelNam, dateStrTime'));
      savePlot(plot_path,filename);
    end 

% RDF plots
    inR = [0; inR];
    density1R = [0; density1R];
    density2R = [0; density2R];
    density3R = [0; density3R];
    density4R = [0; density4R];
    density5R = [0; density5R];
    totalDensityR = [0; totalDensityR];


    figure('Name','RDFcomponents'),
    plot(inR,density1R,'-b','LineWidth',2,'MarkerIndices',1:60:length(density1R));
    hold on;
    plot(inR,density2R,'-r','LineWidth',2,'MarkerIndices',1:60:length(density2R));
    hold on;
    if(model == 5) || (model == 6) || (model == 7)
       plot(inR,density3R,'-c','LineWidth',2,'MarkerIndices',1:60:length(density3R));
       hold on;
    end 

    if(model == 7)
       plot(inR, density5R, 'Color', [1, 0.5, 0], 'LineWidth', 2, 'MarkerIndices', 1:60:length(density5R));
       hold on;
    end 

    plot(inR,density4R,'-g','LineWidth',2,'MarkerIndices',1:60:length(density4R));
    hold on;
    plot(inR,totalDensityR,'-m','LineWidth',2,'MarkerIndices',1:60:length(totalDensityR));
    hold on;
    
    if(model <= 4)
        legend('Short range',...
           '2s',...
           'Long Range',...
           'Total',...
           'Location','NorthEast');
    elseif (model == 5) || (model == 6)
        legend('Short range',...
           '2s',...
           '2p',...
           'Long Range',...
           'Total',...
           'Location','NorthEast');

    elseif (model == 7)
        legend('Short range',...
           '2s',...
           '2p',...
           '3s',...
           'Long Range',...
           'Total',...
           'Location','NorthEast');
    end

    hold off;
    xlabel('Radial Distance, $\rm{r~(au)} $ ','interpreter','LaTex')
    ylabel('$4\pi r^2 \rho(r)$ ','interpreter','LaTex')
    titlename = char(sprintf('Model radial distribution functions'));
    title(titlename);
    filename = char(sprintf('componentRDFs_%s_%s', modelNam, dateStrTime'));
    savePlot(plot_path,filename)

%% Save table of values in output file
% fprintf('Table of parameters used in this fitting routine...\n')
% Create two columns of parameters used for the fitting and number of electrons

end

if(model == 1)
Parameters = ["Atomic state"; "emax"; "A0_dbar"; "Z"; "N"; "alpha"; "B0";...
    "";"Data source"; ...
    "dBar"; "Analytical model"];
Values = [atomicstate; eps; A0_dBar; Z; Ne; alpha; B0;...
    "";"integrated # e-";...
    Int_dBar; Int_model];

T1 = table(Parameters,Values);
head(T1);
writetable(T1,outLog2,'delimiter','\t')
T_dataPlot = [inR, ln_model_fit, radial_model_fit];
writematrix(T_dataPlot, outLog4,'Delimiter','space')

elseif(model == 2)

Parameters = ["Atomic state"; "emax"; "A0_dbar"; "Z"; "N"; "alpha"; "B0";...
    "";"Boundary Conditions"; ...
    "B0LB"; "alphaLB"; "B0UB"; "alphaUB"; ...
    "";"Data source"; ...
    "G09";"dBar"; "Analytical model"];
Values = [atomicstate; eps; A0_dBar; Z; Ne; alpha; B0;...
    "";"BCs";...
    B0LB; alphaLB; B0UB; alphaUB; ...
    "";"integrated # e-";...
    Int_dBar; Int_model];

T1 = table(Parameters,Values);
head(T1);
writetable(T1,outLog2,'delimiter','\t')
T_dataPlot = [ln_model_fit, radial_model_fit];
writematrix(T_dataPlot, outLog4,'Delimiter','space')

elseif(model >= 3)
    
% Note: Z = Zeff in Models 4, 5, and 6.
    if(model == 3)
        Zeff = Z;
        

    end
    
%  save values for imputting as documentation in Excel spreadsheet.
eps_opt   = (alpha^2)/2;               % 1st long-range density parameter
beta_opt  = beta;                      % optimized beta
A0_opt    = A0;
rho0      = model_fit(1,1);            % rho(0) from model_fit after evaluation of model
inA0      = inRho(1,1);                % A0 from data generated either dBar or G09

% ************************************************************************************
KatoNum = Z*A0 + C0*mm*Z + gamma*C0*mm^2 + E0*qq*pp*Z + delta*E0*pp^2;
KatoDen = Z*(A0 + C0*mm^2 + E0*pp^2);
Kato    = KatoNum/KatoDen;
% **********************************************************************************

Parameters = ["Atomic state"; "emax"; "A0_dbar"; "Z"; "N"; "alpha"; "B0";...
    "---";"Boundary Conditions"; ...
    "B0LB"; "alphaLB"; "C0LB"; "gammaLB"; "mmLB"; "betaLB"; "A0LB"; "D0LB"; "etaLB"; 'bbLB'; 'E0LB'; 'deltaLB'; 'ppLB'; 'qqLB'; 'ttLB'; ...
    "B0UB"; "alphaUB"; "C0UB"; "gammaUB"; "mmUB"; "betaUB"; "A0UB"; "D0UB"; "etaUB"; 'bbUB'; 'E0UB'; 'deltaUB'; 'ppUB'; 'qqUB'; 'ttUB';...
    "---";"Input Values"; ...
    "r1"; "r2"; "w1"; "w2"; "w3"; "w4"; "w5"; ...
    "---"; "Initial Values"; ...
    "ini_B0"; "ini_alpha"; "ini_C0"; "ini_mm"; "ini_gamma"; "ini_D0"; "ini_eta"; 'ini_bb'; 'ini_E0'; 'ini_delta'; 'ini_pp'; 'ini_qq'; 'ini_tt';...
    "---"; "Fitted Values"; ...
    "B0"; "alpha"; "C0"; "mm"; "gamma"; "A0"; "beta"; "Zeff"; "D0"; "eta"; 'bb'; 'E0'; 'delta'; 'pp'; 'qq'; 'tt';...
    "---";"Computed values"; ...
    "betaFFAZ"; "Kato"; "Objtv Fn"; 
    "---";"Data source"; ...
    "dBar"; "Analytical model";
    ];

Values = [atomicstate; eps; A0_dBar; Z; Ne; alpha; B0;...
    "---";"BCs";...
    B0LB; alphaLB; C0LB; gammaLB; mmLB; betaLB; A0LB; D0LB; etaLB; bbLB; E0LB; deltaLB; ppLB; qqLB; ttLB;...
    B0UB; alphaUB; C0UB; gammaUB; mmUB; betaUB; A0UB; D0UB; etaUB; bbUB; E0UB; deltaUB; ppUB; qqUB; ttUB;...
    "---";"Ins"; ...
    r1; r2; w1; w2; w3; w4; w5; ...
    "---"; "X0"; ...
    ini_B0; ini_alpha; ini_C0; ini_mm; ini_gamma; ini_D0; ini_eta; ini_bb; ini_E0; ini_delta; ini_pp; ini_qq; ini_tt;...
    "---"; "BestX"; ...
    bestx(1); bestx(2); bestx(3); bestx(4); bestx(5);...
    bestx(6); bestx(7); bestx(8); bestx(9); bestx(10); bestx(11); bestx(12); bestx(13); bestx(14); bestx(15); bestx(16);...
    "---"; "CompVals"; ...
    betaFFAZ; Kato; ObjFn; ...
    "---";"integrated # e-";...
    Int_dBar; Int_model];

T1 = table(Parameters,Values)
head(T1);
writetable(T1,outLog2,'delimiter','\t')

assert(fclose(fid2)==0);  % check for successful file close

% Add new variables at the end of the list, with 99999999 placeholder for Notes
% column in Excel spreadsheet, to maintain back-compatibility as number of
% states (and thus, fitting parameters) is expanded.

NotesPlaceholder = 99999999;

T_excel = [dateNumTime, inZ, Ne, ini_eps, ini_alpha, ini_beta, model, inA0, ini_B0, ini_C0, ini_mm, ini_gamma,...
    mmLB, mmUB, alphaLB, alphaUB, betaLB, betaUB, gammaLB, gammaUB, A0LB, A0UB, B0LB, B0UB, C0LB, C0UB, ...
    r1, r2, w1, w2, w3, w4, w5, ...
    mm, alpha, eps_opt, beta_opt, gamma, A0_opt, B0, C0, Zeff, Int_model, rho0, Kato, ObjFn, ...
    NotesPlaceholder, rCut, ini_D0, D0LB, D0UB, ini_eta, etaLB, etaUB, D0, eta, betaFFAZ, ini_bb, bbLB, bbUB, bb, ...
    ini_E0, E0LB, E0UB, E0, ini_delta, deltaLB, deltaUB, delta, ini_pp, ppLB, ppUB, pp, ini_qq, qqLB, qqUB, qq, ini_tt, ttLB, ttUB, tt];

writematrix(T_excel, outLog3,'Delimiter','space');

% Use padcat to make function vector lengths the same for plotting
T_dataPlot = padcat(ln_model_fit, radial_model_fit, lnrhodBarOver4pi, radial_dBar); 
writematrix(T_dataPlot, outLog4,'Delimiter','space');

end
commandwindow
