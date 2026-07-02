%% Model 2: curve fitting via optimization of B0 and alpha/beta
function ObjFn = Model2ObjFn(x ,r, rhoRef, A0, Z, Ne)

fourpi = 4.0*pi;

B0    = x(1);  
alpha = x(2);  
beta  = ((Z-Ne+1)/alpha - 1);

% Model density is the sum of short-range (F1) and long-range (F2) density functions
F1         = A0.*(exp(-2.*Z.*r));
betav      = (r<=0.2).*1 + (r>0.2).*beta;
F2         = (r>0.2).*(B0.*((r).^(2.*betav)).*(exp(-2.*alpha.*(r))));
den_fit    = F1 + F2;

% Compute integrated density using analytic solution for Model 2 (see intgrl_B0_num
% comments), and current values of alpha, beta, and B0.

I1 = .25*1.0/(Z^3);

fac1 = 2*beta+3;
fac2 = (2.0*alpha)^fac1;
fac2 = 1.0/fac2;
fac3 = gamma(fac1);
I2   = fac2 * fac3;

NeComp = fourpi*(A0*I1 + B0*I2);

% Objective function is the square difference between target and fitted
% model + square difference of exact and model-integrated number of electrons.

% Hardwire weights for now--- bias toward conserving number of electrons
w1 = .05;
w2 = .95;
ObjFn = w1*sum((rhoRef - den_fit).^2) + w2*(Ne - NeComp)^2;
% ObjFn = sum((rhoRef - den_fit).^2);

end


