function Y_vec = MatMulVecOp(A, X_vec)
% performs matrix-vector multiplication between a matrix and a vector
% operator represented by a stack of matrices.
% Operation in LaTex notation:
%                   $$ \hat{\vec{Y} = A \hat{\vec{X}} $$
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% A     : a [p, n] matrix
% X_vec : a [m, m, n] stack of matrices representing the input vector
%         operator. The operators are indexed along the 3rd dimension.
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% Y_vec : a [m, m, p] stack of matrices representing the output vector
%         operator. The operators are indexed along the 3rd dimension.

assert(size(A,2)==size(X_vec,3),'The matrix dimensions are not compatible with dimensions of the vector operator');
Y_vec = squeeze(sum(reshape(A.',[1,1,size(A.')]).*reshape(X_vec,[size(X_vec),1]),3));

end