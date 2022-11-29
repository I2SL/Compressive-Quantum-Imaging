function S_stack = SLD(F_stack, R)
% A function for evaluating S in the implicit matrix equation for the
% Symmetric Logarithmic Derivative
%                       2F = RS + SR
% This method requires that F and R be Hermitian.
%
% Reference: "Quantum Estimation for Quantum Technology" - (Eqn. 12)
% Matteo G. A. Paris (2009)
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% F_stack   : a stack of matrices corresponding to the operator F in the
%             implicit SLD equation
% R         : a matrix corresponding to the operator R in the implicit SLD
%             equation
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% S_stack   : solutions to the implicit SLD equations

% Add dimensionality to F if its not a matrix stack (i.e. just 1 matrix)
is_stack = 1;
if length(size(F_stack)) < 3
    F_stack = reshape(F_stack,[size(F_stack),1]); 
    is_stack = 0;
end


% Check that the inputs are Hermitian
assert(ishermitian(R),'R matrix in implicit SLD equation is not Hermitian');
for j =1:size(F_stack,3)
    assert(ishermitian(F_stack(:,:,j)),['Matrix ',num2str(j),' in F_stack is not Hermitian']);
end

% Get the eigenvectors/values of R
[V_R, d] = eig(R,'vector');

% eigenvalue grids for creating all possible sum pairs
[D1,D2] = meshgrid(d,d);

% Compute an S for each matrix in F_stack
S_stack = 2 * mtimesx(V_R',mtimesx(F_stack,V_R))./(D1+D2 + (D1+D2 == 0) ).*(D1+D2 ~= 0);

% Transform the S matrices back to the original representation space
S_stack = mtimesx(V_R,mtimesx(S_stack,V_R'));

% Average each output the matrix with its Hermitian conjugate for numerical stability
% S_stack = 1/2 * (S_stack + mtimesx(S_stack,'C',eye(size(S_stack,1,2))));

% check that the output is Hermitian
for j =1:size(S_stack,3)
    % Average each output the matrix with its Hermitian conjugate for numerical stability
    S_stack(:,:,j) = 1/2 * (S_stack(:,:,j) + S_stack(:,:,j)'); 
    assert(ishermitian(S_stack(:,:,j)),['Matrix ',num2str(j),' in S_stack is not Hermitian']);
end

% Remove added dimension if F was not a stack
if ~is_stack
    S_stack = squeeze(S_stack);
end

end