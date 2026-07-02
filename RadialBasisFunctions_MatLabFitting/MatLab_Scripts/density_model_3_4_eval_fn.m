%% Model fit to average density, using expression of Wang and Parr. See dissertation text for details.
%  Note: Fewer ranges needed for Period 2 elements (K, L shells only)
%        May not want to use ranges as in Wang and Parr-- just 2s-like
%        functional form to capture shell structure (added to Model 2 short-range
%        and long-range forms).
%
%  Godwin Amo-Kwao, 7/25/19.
%

function den_fn = density_model_3_4_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, r, rCut)

% beta   = ((Z-Ne+1)/alpha - 1);

F1      = A0*(exp(-2*Z*r));                                     % short-range 
betav   = (r==0).*1 + (r>0).*beta;

% -------------------------------------------------------------------
F2_a   = (mm  - Z*r);
F2_b   = F2_a .* F2_a;
F2     = C0*F2_b.*exp(-2*gamma*r);                              % 2s shell structure        

% set rCut = 0.0 to default to previous behavior.
F3     = (r>rCut).*(B0*(r.^(2*betav)).*(exp(-2*alpha*r)));      % long-range

% Model density is the sum of short-range (F1), long-range (F3), and 2s-shell (F2) density functions
den_fn = F1 + F2 + F3;

end