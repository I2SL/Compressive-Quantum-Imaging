classdef MVNPrior < Prior
    properties
        % mean and covariance
        mu
        cov
        
        % compositional flag
        isCompositional = 0
        
    end
    
    methods
        % constructor
        function obj = MVNPrior(params)
            % params should be a struct with fields:
            obj.mu = params.mu;
            obj.cov = params.cov;
        end
        
        % prior probability
        function [p, ln_p, grad_ln_p] = prob(obj,x)
            
            % prior
            p = mvnpdf(x',obj.mu',obj.cov)';

            % log prior
            xx = x - obj.mu; % center samples by mean
            ln_p = - 0.5 * sum(xx .* (obj.cov \ xx), 1); % + const.

            % gradient of log prior
            grad_ln_p = - obj.cov \ xx;
            
        end
        
        % random sampler
        function x = rnd(obj,n)
            x = mvnrnd(obj.mu',obj.cov,n)';
        end
        
        % match first and second moments
        function obj = match(obj,x_mu,x_cov) 
            obj.mu = x_mu;            
            obj.cov = x_cov;
        end
    end
        
    
end
