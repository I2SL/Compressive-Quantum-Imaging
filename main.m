function main(array_id)
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% array_id          - array id for cluster computing and random number generation

seed = str2double(array_id);

n_post_samples = 2e5;       % number of posterior samples
n_HG_samples = (2^8)+1;     % number of samples of the HG polynomials
n_HG_modes = 15;             % number of HG modes
n_photons =  5e3;           % number of photons collected per adaptive measurement
max_iter = 50;              % max number of bayesian iteration
depth = 2;                  % depth of wavelet decomposition
rl = 0.2;                   % rayleigh length (in fractions of signal support)

[MMSE_est, mu_seq, cov_seq, l_seq, B_seq, x0, T, W] = simulateCQE(n_post_samples, n_photons, max_iter,...
                                                                'rl',rl,...
                                                                'n_HG_samples',n_HG_samples,...
                                                                'n_HG_modes',n_HG_modes,...
                                                                'prior_name','gbm',...
                                                                'posterior_method', 'conjugate',...
                                                                'rng_seed',seed,...
                                                                'depth',depth....
                                                                );

save

end