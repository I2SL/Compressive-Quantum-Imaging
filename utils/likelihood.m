function [p_like, ln_like, grad_ln_like] = likelihood(x, l_vec_set, B_gamma_set, W_vec)
% The likelihood probability of observing l_vec under the measurement
% operator B_gamma (the joint estimator). The likelihood function is a
% multinomial over the probabilty distribution of the outcomes.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% l_vec_set   : a matrix whose columns are measurement vectors containing 
%               the number of photons detected in each of the eigenstates 
%               of the measurement operators in B_gamma_set
% B_gamma_set : a stack matrices representing a sequence of measurement operators
% x           : a matrix whose columns are samples of the parameter vector
% W_vec       : the parameter operators 
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% p_l         : the cumulative probability of the outcomes in l_vec_set
% ln_like     : the log of the cumulative probability (up to an additive constant)
% grad_ln_like: the gradient of the log probability

% initialize outputs
p_like = 1;
ln_like = 0;
grad_ln_like = 0;

% loop through each measurement and compute its contribution to the cumulative likelihood
for i = 1:size(B_gamma_set,3)
    [p_i,ln_i,grad_ln_i] = likelihood_single(x,l_vec_set(:,i),B_gamma_set(:,:,i),W_vec);
    p_like = p_like .* p_i;
    ln_like = ln_like + ln_i;
    grad_ln_like = grad_ln_like + grad_ln_i;
end
    
end


function [p_like, ln_like, grad_ln_like] = likelihood_single(x, l_vec, B_gamma, W_vec)
% The likelihood probability of observing l_vec under the measurement
% operator B_gamma (the joint estimator). The likelihood function is a
% multinomial over the probabilty distribution of the outcomes.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% l_vec       : measurement vector containing the number of photons detected in each
%               of the eigenstates of the joint estimator B_gamma l_vec = [n_1,n_2,...,n_N]
% B_gamma     : the measurement matrix for the joint estimator.
% x           : a matrix whose columns are samples of the parameter vector
% W_vec       : the parameter operators 
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% p_l         : the probability of the outcome l_vec
% ln_like     : the log probability of the outcome l_vec (up to an additive constant)
% grad_ln_like: the gradient of the log probability


% dimension constants
N = size(x,1);       % number of parameters
M = size(x,2);       % number of parameter instances for which to compute likelihood
D = size(B_gamma,1);  % dimensionality of hilbert space

% get eigenstates of measurement operator
[V,~] = eig(B_gamma);

%% compute the likelihood
% container for the eigenstate probabilities
p_eig = zeros([D+1,M]);

% super fast way of computing eigenstate probabilities for each instantiation 
% of the sampled density operators (rhos) without breaking memory limits.
WW = squeeze(sum(mtimesx(V',W_vec).*V.',2));
p_eig(1:D,:) = WW*x; 


% throw a warning if some of the probabilities are negative 
% (suggests a numerical instability in eigen decomposition of operators)
neg_prob = p_eig<0;
if any(neg_prob(:))
    warning('Negative probabilties computed for modal measurement outcomes in likelihood calculation. Setting negative probabilities to zero. Consider up-sampling the image plane.')
    p_eig(neg_prob) = 0;
end

% compute the residual probability of photon counts in residual mode
residual = 1-sum(p_eig,1);
residual(residual < 0) = 0;

% insert the residual probability
p_eig(D+1,:) = residual;

% normalize probabilities
p_eig = p_eig./sum(p_eig,1);

% likelihood is a multinomial
p_like = mnpdf(l_vec', p_eig')';

%% compute log likelihood
rep_l_vec = repmat(l_vec,[1,M]);
ln_like = sum(rep_l_vec.*real(log(p_eig + (rep_l_vec == 0 & p_eig == 0))),1); % + const.

%% compute gradient of log likelihood
% compute coefficients for operator expectations
lp_ratio = l_vec./p_eig;


% zero-out elements where photon counts or eigenstate probabilities were 0
rmv = (isnan(lp_ratio) | isinf(lp_ratio));
lp_ratio(padarray(rmv(1:D,:),[1,0],0,'post')) = 0;

% compute global coefficient prefactor
prefactor = zeros(1,M);
prefactor(~rmv(D+1,:)) = 1-lp_ratio(D+1,~rmv(D+1,:));

% gradient of log likelihood
grad_ln_like = prefactor .* (WW.'*lp_ratio(1:D,:));

end