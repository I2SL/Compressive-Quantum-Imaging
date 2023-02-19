function [MMSE_est,...
          x_mu_seq,  x_cov_seq, ...
          l_vec_seq, B_gamma_seq] = AdaptiveMeasurement(B_start, Y_vec, posterior_obj, measurement_obj)
% Executes a bayesian inference protocol with adaptive measurements
% for estimating the parameters of a quantum density state.
%
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% B_start           : the initial measurement operator
% Y_vec             : a stack of parameter operators
% posterior_obj     : an instance of a posterior object
% measurement_obj   : an instance of a measurement object
% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% MMSE_est          : the Minimum-Mean-Squared Error estimate of the
%                     parameters

% iteration
i = 1;

% initialize parameter moment sequences
x_mu_seq = [];
x_cov_seq = [];

% initialize measurement sequences
[measurement_obj,l_start] = measurement_obj.MakeMeasurement(B_start);
l_vec_seq(:,i) = l_start;
B_gamma_seq(:,:,i) = B_start;


% run adaptive Bayesian measurement scheme
while ~measurement_obj.stop_flag     
    
    % display iteration
    disp(i);

    % update the parameter moments given new measurement
    [posterior_obj, x_mu,x_cov] = posterior_obj.PosteriorMoments(l_vec_seq, B_gamma_seq, Y_vec,...
                                                  x_mu_seq, x_cov_seq);
                                              
    % compute the next measurement operator                                        
    if posterior_obj.prior_obj.isCompositional
       B_gamma = JointParameterEstimator(Y_vec, x_mu, x_cov);
    else
       B_gamma = JointParameterEstimator(Y_vec(:,:,2:end),x_mu(2:end),x_cov(2:end,2:end)); % crop out augmented parameter
    end

    % get the next measurement 
    [measurement_obj,l_vec] = measurement_obj.MakeMeasurement(B_gamma);

    % increment iteration count
    i = i + 1;

    % add adaptive results to sequences
    x_mu_seq(:,i) = x_mu;
    x_cov_seq(:,:,i) = x_cov;
    l_vec_seq(:,i) = l_vec;
    B_gamma_seq(:,:,i) = B_gamma;
    
    save('x_mu_seq','x_cov_seq','l_vec_seq','B_gamma_seq')

end

% return parameter Minimum Mean-Squared Error Estimator 
MMSE_est = x_mu; 

end



