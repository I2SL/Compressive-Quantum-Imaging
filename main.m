function main(array_id)
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% array_id          - array id for cluster computing and random number generation

seed = str2double(array_id);

n_post_samples = 1e5;
n_photons =  1e5;
max_iter = 50;

[MMSE_est, mu_seq, cov_seq, l_seq, B_seq, x0, T, W] = simulateCQE(n_post_samples, n_photons, max_iter,'rng_seed',seed);

save

%{

% show target params vs reconstruction
figure
x0 = (M*V)\w0;
stem(x0,'filled','black')
hold on
stem(MMSE_est,'filled','red')
hold off

%}

end