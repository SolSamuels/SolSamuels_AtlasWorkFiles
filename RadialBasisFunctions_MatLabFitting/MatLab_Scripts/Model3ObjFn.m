%% Model 3: curve fitting via optimization of B0, alpha/beta, C0, mm, and gammaf

function ObjFn = Model3ObjFn(x, r, rhoRef, Z, Ne, w1, w2, w3, w4, w5, r1, r2, rCut)

fourpi = 4.0*pi;

B0     = x(1);  
alpha  = x(2);  
C0     = x(3);
mm     = x(4);
gammaf = x(5);
A0     = x(6);
beta   = x(7);

% Model density is the sum of short-range (F1), long-range (F3), and 2s-shell (F2) density functions
F1     = A0*(exp(-2*Z.*r));

%----------------------------------------------------------------------

F2_a       = mm-(Z*r);
F2_b       = (F2_a).*(F2_a);
F2         = C0*F2_b.*exp(-2*gammaf*r);

% set rCut = 0.0 to reproduce previous behavior.
betav  = (r==0).*1 + (r>0).*beta;
F3         = (r>rCut).*(B0*((r).^(2.*betav)).*(exp(-2*alpha*r)));

%---------------------------------------------------------------------------
% Model density is the sum of short-range (F1), long-range (F3), and 2s-shell (F2) density functions
den_fit    = F1 + F2 + F3;

% Compute integrated density (4*pi*r^2*den_fit) for applying total # electron constraint

%----First term, short-range----
fun1 = @(r)r.^2.*exp(-2*Z*r);
int_fun1 = integral(fun1,0,Inf);

% check: compare numerical integration result int_fun1 with analytic result I1
I1 = .25*1.0/(Z^3);

%----Second term, long-range----
fun2 = @(r)(r.^(2*beta + 2).*exp(-2*alpha*r));
int_fun2 = integral(fun2,0,Inf);

% check: compare numerical integration result int_fun2 with analytic result I2
fac1 = 2*beta+3;
fac2 = (2.0*alpha)^fac1;
fac2 = 1.0/fac2;
fac3 = gamma(fac1);
I2   = fac2 * fac3;

%----Third term, 2s shell structure----
fun3     = @(r)(r.^2).*((mm-(Z*r)).^2).*exp(-2*gammaf*r);
int_fun3 = integral(fun3,0,Inf);

% check: compare numerical integration result int_fun3 with analytic result I3

gammacub   = gammaf^3;
gammaquar  = gammacub*gammaf;
gammaquint = gammaquar*gammaf;

J1 = (mm*mm)/(4*gammacub);
J2 = -(.75)*Z*mm/gammaquar;
J3 = .75*Z^2/gammaquint;

I3 = J1 + J2 + J3;

%% Function expression for r^2*den(r)
rsq_den_fit=@(r)r.^2.*(A0.*(exp(-2.*Z.*r)) +...
    (r>0).*(B0.*((r).^(2.*((r==0).*1 + (r>0).*beta))).*(exp(-2.*alpha.*(r)))) +...
    C0.*((mm-(Z*r)).^2).*exp(-2*gammaf*r));

NeComp = fourpi*integral(rsq_den_fit,0,Inf);

% check: compare numerical integration result Necomp with analytic result Necomp_check 
Necomp_check = fourpi*(A0*I1 + B0*I2 + C0*I3);

% Objective function is the weighted sum of square difference between target and fitted
% model, square difference of exact and model-integrated number of electrons, and 
% square difference of exact and model-computed Kato ratio.

radialRef       = fourpi*r.^2.*rhoRef;           % RDF for reference data
radialFit       = fourpi*r.^2.*den_fit;          % RDF for analytic model with fitted B0 

DiffSq          = (radialRef - radialFit).^2;
mask1           = (r <= r1);
mask2           = (r > r1 & r <= r2);
mask3           = (r > r2);

N1 = sum(mask1); 
N2 = sum(mask2);
N3 = sum(mask3); 

msd1 = sum(DiffSq.*mask1)/N1;
msd2 = sum(DiffSq.*mask2)/N2;
msd3 = sum(DiffSq.*mask3)/N3;

KatoNum = Z*A0 + C0*mm*Z + gammaf*C0*mm^2;
KatoDen = Z*(A0 + C0*mm^2);
Kato    = KatoNum/KatoDen;

ObjFn =  w1*msd1 + w2*msd2 + w3*msd3 + w4*(Ne - NeComp)^2 + w5*(1 - Kato)^2;

end
