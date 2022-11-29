function T0 = Gamma_0(Y_vec,x_mu)
% Computes the mean density operator.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% Y_vec    : a stack of matrices representing the parameter operators
% x_mu     : expectation of the parameters
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% T0       : the mean density operator

T0 = MatMulVecOp(x_mu.',Y_vec);
    
end