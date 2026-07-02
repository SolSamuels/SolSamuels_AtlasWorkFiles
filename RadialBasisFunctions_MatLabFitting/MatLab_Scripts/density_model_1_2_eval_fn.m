%% Model fit to average density, expressed as the sum of a short-range (Kato cusp exponential) term
% and long-range asymptotic term, where the exponential decay depends on the IP of the atom or ion 
% whose density is being modeled.  See dissertation text for details.
%
%  Godwin Amo-Kwao, 7/12/19.
%
%  Modified to deal with long-range term blowup at r=0 for anions.
%
%  Version 1.5, SRA, 7/26/19
%
%  Original if statement syntax did not work:
%    if r <= 0.01
%       den_fn = A0*(exp(-2*Z*r));
%    else
%       den_fn = A0.*(exp(-2.*Z.*r)) + B0_integral_num.*((r).^(2.*beta))...
%         .*(exp(-2.*alpha.*r));
%
%  Re-expressed this way:
%     F1  = A0.*(exp(-2.*Z.*r))
%     F2a = (r>0).*( B0.*((r).^(2.*beta)))
%     F2b = (r>0).*(exp(-2.*alpha.*(r)))
%     F2  = F2a.*F2b
%     den_fn = F1+F2

function den_fn = density_model_1_2_eval_fn(A0, B0, Z, alpha, Ne, r)

beta  = ((Z-Ne+1)/alpha - 1);

F1 = A0.*(exp(-2.*Z.*r));               % short-range 
beta = (r==0).*1 + (r>0).*beta;


F2 = (r>0).*(B0.*((r).^(2.*beta)).*(exp(-2.*alpha.*(r))));  % long-range

% Model density is the sum of short-range (F1), long-range (F2),
den_fn = F1 + F2;

end
