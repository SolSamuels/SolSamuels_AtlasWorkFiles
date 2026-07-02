%% Model 6: curve fitting via optimization of B0, alpha/beta, C0, mm, gammaf, A0, Zeff, D0, eta, and bb

function ObjFn = Model6ObjFn(x, r, rhoRef, Ne, w1, w2, w3, w4, w5, r1, r2, rCut)

fourpi = 4.0*pi;

B0     = x(1);  
alpha  = x(2);  
C0     = x(3);
mm     = x(4);
gammaf = x(5);
A0     = x(6);
beta   = x(7);
Zeff   = x(8);
D0     = x(9);
eta    = x(10);
bb     = x(11);

% Model density is the sum of short-range (F1), 2s-shell (F2), 2p-shell (F3), and 
% long-range (F4) density functions

%---------------------------------------------------------------------------

% Compute density components
F1 = A0*exp(-2*Zeff.*r);
F2 = C0*(mm-(Zeff*r)).^2 .* exp(-2*gammaf*r);
F3 = D0*(Zeff*r).^2 .* exp(-2*eta*r);

betav = (r==0).*1 + (r>0).*beta;
betav = betav + bb*r;  % linear beta
F4 = (r>rCut) .* (B0*(r.^(2.*betav)) .* exp(-2*alpha*r));

den_fit = F1 + F2 + F3 + F4;

% Clip to prevent NaNs and Infs
den_fit(~isfinite(den_fit)) = 0;
den_fit = max(den_fit, 1e-12);
den_fit = min(den_fit, 1e12);

% Compute integrated density (4*pi*r^2*den_fit) for applying total # electron constraint

%----First term, short-range----
fun1 = @(r)r.^2.*exp(-2*Zeff*r);
int_fun1 = integral(fun1,0,Inf);

% check: compare numerical integration result int_fun1 with analytic result I1
I1 = .25*1.0/(Zeff^3);

% Comment-out print statements for production calculations
% fprintf(1,'Numerical I1= %f.\n',int_fun1);
% fprintf(1,'Analytic  I1= %f.\n',I1);

%----Second term, long-range----

% Original:  fun2 = @(r)(r.^(2*beta + 2).*exp(-2*alpha*r));
% "Quad"     fun2 = @(r)(r.^(2*((r==0).*1 + (r>0).*beta + bb.*r + qq.*r.*r) + 2).*exp(-2*alpha.*r));
fun2 = @(r)(r.^(2*((r==0).*1 + (r>0).*beta  + bb.*r) + 2).*exp(-2*alpha*r));
int_fun2 = integral(fun2,0,Inf);

% check: compare numerical integration result int_fun2 with analytic result I2
% Note: this is no longer correct for r-dependent \beta in Model 6--- need
% to fix.
% fac1 = 2*beta+3;
% fac2 = (2.0*alpha)^fac1;
% fac2 = 1.0/fac2;
% fac3 = gamma(fac1);
% I2   = fac2 * fac3;

% Comment-out for production calculations
% fprintf(1,'Numerical I2 (LR)= %f.\n',int_fun2);
% fprintf(1,'Analytic  I2 (LR)= %f.\n',I2);

%----Third term, 2s shell structure----
fun3     = @(r)(r.^2).*((mm-(Zeff*r)).^2).*exp(-2*gammaf*r);
int_fun3 = integral(fun3,0,Inf);

% check: compare numerical integration result int_fun3 with analytic result I3

gammacub   = gammaf^3;
gammaquar  = gammacub*gammaf;
gammaquint = gammaquar*gammaf;

J1 = (mm*mm)/(4*gammacub);
J2 = -(.75)*Zeff*mm/gammaquar;
J3 = .75*Zeff^2/gammaquint;

I3 = J1 + J2 + J3;

% Comment-out for production calculations
% fprintf(1,'Numerical I3 (2s) = %f.\n',int_fun3);
% fprintf(1,'Analytic  I3 (2s) = %f.\n',I3);

%----Fourth term, 2p shell structure----
fun4     = @(r)(r.^2).*((Zeff*r).^2).*exp(-2*eta*r);
int_fun4 = integral(fun4,0,Inf);

etaquint = eta^5;
I4       = .75*Zeff^2/etaquint;

% Comment-out for production calculations
% fprintf(1,'Numerical I4 (2p) = %f.\n',int_fun4);
% fprintf(1,'Analytic  I4 (2p) = %f.\n',I4);

%% Function expression for r^2*den(r)
% rsq_den_fit=@(r)r.^2.*(A0.*(exp(-2.*Zeff.*r)) +...
%     (r>0).*(B0.*((r).^(2.*((r==0).*1 + (r>0).*beta))).*(exp(-2.*alpha.*(r)))) +...
%     C0.*((mm-(Zeff*r)).^2).*exp(-2*gammaf*r) + D0.*(Zeff*r).^2.*exp(-2*eta*r));

rsq_den_fit=@(r)r.^2.*(A0.*(exp(-2.*Zeff.*r)) +...
    (r>0).*(B0.*((r).^(2.*((r==0).*1 + (r>0).*(beta + bb.*r)))).*(exp(-2.*alpha.*(r)))) +...
    C0.*((mm-(Zeff*r)).^2).*exp(-2*gammaf*r) + D0.*(Zeff*r).^2.*exp(-2*eta*r));

NeComp = fourpi*integral(rsq_den_fit,0,Inf);

% check: compare numerical integration result Necomp with analytic result Necomp_check 
% NeComp_check = fourpi*(A0*I1 + B0*I2 + C0*I3 + D0*I4);

% fprintf(1,'Numerical Ne = %f.\n',NeComp);
% fprintf(1,'Analytic  Ne = %f.\n',NeComp_check);

% Objective function is the weighted sum of square difference between target and fitted
% model, square difference of exact and model-integrated number of electrons, and 
% square difference of exact and model-computed Kato ratio.

% Compute radial distribution functions
radialRef       = fourpi*r.^2.*rhoRef;           % RDF for reference data
radialFit       = fourpi*r.^2.*den_fit;          % RDF for analytic model with fitted B0 

% Compute objective MSD in three regions
DiffSq = (radialRef - radialFit).^2;
mask1 = (r <= r1);
mask2 = (r > r1 & r <= r2);
mask3 = (r > r2);

N1 = max(1,sum(mask1));
N2 = max(1,sum(mask2));
N3 = max(1,sum(mask3));

msd1 = sum(DiffSq.*mask1)/N1;
msd2 = sum(DiffSq.*mask2)/N2;
msd3 = sum(DiffSq.*mask3)/N3;

% Integrated electron count using trapz 
%NeComp = trapz(r, radialFit);

% Kato condition with safe denominator
KatoNum = Zeff*A0 + C0*mm*Zeff + gammaf*C0*mm^2;
KatoDen = max(Zeff*(A0 + C0*mm^2), 1e-12);
Kato    = KatoNum/KatoDen;

% Weighted objective
ObjFn = w1*msd1 + w2*msd2 + w3*msd3 + w4*(Ne - NeComp)^2 + w5*(1 - Kato)^2;

% Final safety clamp
if ~isfinite(ObjFn)
    ObjFn = 1e12;
end

end
