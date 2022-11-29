function [MMSE_est, mu_seq, cov_seq, l_seq, B_seq, x0, T, W] ...
         = simulateCQE(n_post_samples, n_photons, max_iter, varargin)
% Simulates Compressive Quantum Estimation for a provided target signal.
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS 
% -------------------------------------------------------------------------
% n_post_samples    - number of posterior samples
% n_photons         - number of photons in each measurement
% max_iter          - max number of bayesian iterations
%
% -------------------------------------------------------------------------
% OPTIONAL INPUTS 
% -------------------------------------------------------------------------
% rl                - Rayleigh length (fractional size of the signal support)
% n_HG_modes        - number of HG modes (along 1D)
% n_HG_samples      - number of samples to discretize HG modes (along 1D)
% prior_name        - name of sparsity prior to apply to parameters [dirichlet,mvn,gbm,laplace,cifar10gbm]
% posterior_method  - name of posterior method to use [K_informative, cumulative, conjugate]
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% MMSE_est          -
% mu_seq            -
%
% 
%
% Author: Nico Deshler
% -------------------------------------------------------------------------


% add subdirectories to search path
addpath('utils')
addpath('classes')
addpath('Adaptive-Bayesian-Quantum-Estimation')

% parse inputs
p = inputParser;

% required parameters
addRequired(p,'n_post_samples');
addRequired(p,'n_photons');
addRequired(p,'max_iter');

% optional parameters
addOptional(p,'target',[]);
addOptional(p,'rl',1);
addOptional(p,'n_HG_modes',12);
addOptional(p,'n_HG_samples',(2^5)+1);
addOptional(p,'prior_name','cifar10gbm');
addOptional(p,'posterior_method','K_informative');
addOptional(p,'K',5);
addOptional(p,'rng_seed',0);
addOptional(p,'d',2);                           % dimensionality of the signal
addOptional(p,'depth',2);                       % wavelet decomposition depth

parse(p, n_post_samples, n_photons, max_iter, varargin{:})


% set the random number generator
rng(p.Results.rng_seed)

% instantiate prior object
prior_obj = setup_prior(p.Results.prior_name);

% instantiate reference distribution
ref_obj = prior_obj;

% instantiate posterior object
post_obj = Posterior(prior_obj,ref_obj);
post_obj.num_samples = n_post_samples;
post_obj.method = p.Results.posterior_method;
post_obj.K = p.Results.K;

% get wavelet matrix
d = p.Results.d;
W = db1WaveletMatrix(p.Results.depth, d);

% get wavelet operators
rl = p.Results.rl;
n_HG_modes = p.Results.n_HG_modes;
n_HG_samples = p.Results.n_HG_samples;
W_vec = W_vec_HG(W,rl,n_HG_modes,n_HG_samples,d);

% energy of each basis function
trW = cell2mat(arrayfun(@(i)trace(W_vec(:,:,i)), 1:size(W_vec,3),'UniformOutput',false)'); % applies trace to each layer in matrix stack
trW(sum(W,1)==0) = 0;  % match zeros of wavelet energy

% constraint matrices
M = M_matrix(trW);          % normalization constraint matrix
V = V_matrix(W,M);          % non-negativity constraint matrix

% set compound constraint matrix based on the prior
if prior_obj.isCompositional
    T = M*V;    % if parameters are compositional then both the trace(rho)=1 and target >= 0 constraints are enforced
else
    T = M;      % otherwise only the trace(rho)=1 constraint is enforced    
end

% transformed wavelet operators
Y_vec = MatMulVecOp(T.',W_vec);

% randomly sample the target from the prior if none is provided
target = p.Results.target;
if isempty(target)
    x0 = prior_obj.rnd(1);                  % sample the prior to generate a target signal
    if ~prior_obj.isCompositional
        % shift non-compositional samples back into constraint polytope
        temp = V\x0; temp = temp - min(temp);   
        x0 = V*temp/sum(temp);
    end
    target = W*T*x0;
    dims = (numel(target)^(1/d)) * ones(1,d);
    target = reshape(target, dims);
end

% make sure target signal is valid
assert(abs(sum(target(:))- 1) < 1e-15, 'target intensity signal not normalized');
assert(all(target(:) >= -1e-15), 'target intensity signal has negative values');

% instantiate measurement object
x0 = (W*T) \ target(:);       % ground truth parameters
rho = MatMulVecOp(x0',Y_vec);
meas_obj = SimulatedMeasurement(rho, n_photons, max_iter);

% instantiate starting measurement from direct detection estimate of target 
dd_target = imgaussfilt(target,rl*size(target));    % filter the target with gaussian PSF to emulate direct detection
x_start = (W*T) \ dd_target(:);                     % set the first modal measurement from the MLE wavelet coeffs recovered from direct detection
x_cov_start = prior_obj.cov;                     % and the covariance matrix of the prior
B_start = JointParameterEstimator(Y_vec, x_start, x_cov_start); % joint parameter estimator from direct detection

% run the adaptive modal measurement protocol
[MMSE_est,mu_seq,cov_seq,l_seq,B_seq] = AdaptiveMeasurement(B_start, Y_vec, post_obj, meas_obj);

end