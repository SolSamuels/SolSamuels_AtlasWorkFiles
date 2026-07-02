%  Model 5 fit to spherically-averaged density.  See GAK dissertation text for details.
%  
%  Baseline version, Godwin Amo-Kwao, 7/25/19.
%  Modified to add 2p radial density to complete Period 2 model and accommodate 
%  negative ions, SR Atlas, 10/5/2020


function den_fn = density_model_5_eval_fn(A0, B0, Z, alpha, beta, C0, mm, gamma, D0, eta, r, rCut)

% ---------------------------------------------------------------------------------------

F1      = A0*(exp(-2*Z*r));                                     % short-range cusp density

F2_a   = (mm  - Z*r);
F2_b   = F2_a .* F2_a;
F2     = C0*F2_b.*exp(-2*gamma*r);                              % 2s shell structure   

F3_a   = Z*r;
F3_b   = F3_a .* F3_a;
F3     = D0*F3_b.*exp(-2*eta*r);                                % 2p shell structure    

% set rCut = 0.0 to default to previous behavior.
betav  = (r==0).*1 + (r>0).*beta;
F4     = (r>rCut).*(B0*(r.^(2*betav)).*(exp(-2*alpha*r)));      % long-range asymptotic density

% ---------------------------------------------------------------------------------------

% Model density is the sum of short-range (F1), 2s-shell (F2), 2p-shell (F3), and  
% long-range (F4) density functions

den_fn = F1 + F2 + F3 + F4;

end