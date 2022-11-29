function E_Q = Sigma_Q(T0, B_vec, x_mu, x_cov)
% Computes the Bayesian Quantum Cramer-Rao Lower Bound matrix
% which is similar to the covariance matrix of the estimators.
%
% Reference: "Bayesian multiparameter Quatnum Metrology with Limited Data"
% Rubio and Dunningahm (2020)
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% T0        : the expectation of the density operator
% B_vec     : a matrix stack of the Personick estimators for each parameter
% x_mu      : expectation of the parameters
% x_cov     : covariance matrix of the parameters
% 
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% E_Q         : the BQCRLB

% get the second-moment matrix
E_xx = x_cov + x_mu*x_mu'; 

g = size(B_vec,3);
K = zeros([g,g]);

for i = 1:g
    for j = 1:g
    K(i,j) = trace(T0 * (B_vec(:,:,i)*B_vec(:,:,j)+B_vec(:,:,j)*B_vec(:,:,i))/2);
    end
end

% Quantum Bayesian Cramer-Rao Lower Bound (QBCRLB)
E_Q = E_xx - K;
end