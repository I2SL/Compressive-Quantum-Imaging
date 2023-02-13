function [B_gamma, E_Q] = JointParameterEstimator(Y_vec, x_mu, x_cov)
% calculates the Personick minimum-mean-squared-error estimator B_gamma 
% for the joint-parameter and the Quantum Bayesian Cramer-Rao Lower Bound 
% (QBCRLB) matrix E_Q.
%
% Reference: "Quantum-inspired Multi-Parameter Adaptive Bayesian Estimation for Sensing and Imaging"
% Kit Lee et. al. (2021)
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% Y_vec    :  a stack of matrices representing the parameter operators
% x_mu     :  expectation of the parameters
% x_cov    :  covariance matrix of the parameters
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% B_gamma  : a hermitian matrix representing the MMSE joint parameter
%            estimator (measurement that collectively reduces total
%            uncertainty in the parameters in the Bayesian setting)
% E_Q      : the Quantum Bayesian Cramer-Rao Lower Bound matrix

% compute Gamma_0 and {Gamma_i1}'s    
T0 = Gamma_0(Y_vec,x_mu);
T1_vec = Gamma_1_vec(Y_vec,x_mu,x_cov);

% compute the Personick MMSE estimators for each parameter
B_vec = SLD(T1_vec, T0);

% compute the QBCRLB matrix
E_Q = Sigma_Q(T0,B_vec,x_mu,x_cov);

% get the eigenvectors of QBCRLB matrix
[V_Q,D_Q] = eig(E_Q);

% get the index of the minimum eigenvector
[~,i] = min(diag(D_Q));

% choose joint-parameter projection direction to be the min eigenvector
h = V_Q(:,i);

% update the joint parameter estimator (measurement matrix) for the subsequent measurement    
B_gamma = MatMulVecOp(h',B_vec);

end
