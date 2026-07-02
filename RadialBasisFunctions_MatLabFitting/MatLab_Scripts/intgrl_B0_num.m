%% This function calculates B0 used in the model density_model_eval_fn

function B0 = intgrl_B0_num(A0, Z, alpha, N)
%  supply A0, Z, alpha, beta, N from densityFitting script.

beta  = ((Z-N+1)/alpha - 1);

fun1 = @(r)r.^2.*exp(-2*Z*r);
int_fun1 = integral(fun1,0,Inf);

% check: compare numerical integration result int_fun1 with analytic result I1
I1 = .25*1.0/(Z^3);
 
fun2 = @(r)r.^(2*beta + 2).*exp(-2*alpha.*r);
int_fun2 = integral(fun2,0,Inf);

% check: compare numerical integration result int_fun2 with analytic result I2
fac1 = 2*beta+3;
fac2 = (2.0*alpha)^fac1;
fac2 = 1.0/fac2;
fac3 = gamma(fac1);
I2 = fac2 * fac3;

B0 = (N - 4*pi*A0*int_fun1)/(4*pi*int_fun2);
 
end
