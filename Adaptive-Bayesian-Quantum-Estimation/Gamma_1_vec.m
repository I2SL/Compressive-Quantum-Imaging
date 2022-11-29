function T1_vec = Gamma_1_vec(Y_vec,x_mu,x_cov)
% Computes the stack of first-moment Personick operators
% for each the parameters.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% Y_vec    :  a stack of matrices representing the parameter operators
% x_mu     :  expectation of the parameters
% x_cov    :  covariance matrix of the parameters
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% T1_vec   : a matrix stack stack of first-moment Personick operators
%            for each parameter --> T1_i = T1_Vec(:,:,i) 
    
% The second moment matrix E[x * x']
E_xx = x_cov + x_mu*x_mu';

% Transform the parameter operator vector by the second moment matrix
T1_vec = MatMulVecOp(E_xx, Y_vec);

end


