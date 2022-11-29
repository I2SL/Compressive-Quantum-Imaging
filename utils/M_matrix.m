function M = M_matrix(tr_A) 
% Computes the matrix for transforming the unconstrained parameters 
% into the wavelet coefficients. The matrix is
% defined such that it enforces the hyperplane constraint imposed by the 
% trace-1 property of the density operator.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% tr_A      : A column vector containing the trace of each wavelet operator 
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% M         : the transform matrix


% get hyperplane constraint matrix (constraint 1: normalized intensity)
ff = tr_A/(tr_A'*tr_A);                     % normalize by norm squared


N = numel(ff);      % Number of parameters
idx = 1:N;          % Index list
S0 = idx(ff==0);    % Index list of 0-valued f' elements
S1 = idx(ff~=0);    % Index list of non-0-valued f' elements
M = zeros([N,N]);   % Instantiate the hyperplane constraint matrix (theta = M * aa)
M(:,1) = ff;        % Set the first column to f'  


% M matrix columns for non-0-valued f' elements
for k = 1:numel(S1)-1
    M(S1(k),k+1) = ff(S1(k+1));
    M(S1(k+1),k+1) = -ff(S1(k));
end

% W matrix columns for 0-valued f' elements
for i = 1:numel(S0)
    j = numel(S1) + i;
    M(S0(i),j) = 1;
end


% Check M matrix conditions
zero_tol = 1e-15;       % tolerance for numerical stability
assert(isequal(M(:,2:end)'*ff < zero_tol,ones(N-1,1)),'columns of M are not orthogonal f vector') % columns 2:N of M are orthogonal to ff
assert(det(M) ~= 0,'M is not invertible')                                                         % M is invertible

end