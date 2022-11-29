classdef DeltaPrior < Prior
    properties
        
        % mean and covariance
        mu
        cov
        
        % compositional flag
        isCompositional
    end
    
    methods
        function obj = DeltaPrior(params)
            obj.mu = params.x;
            obj.cov = zeros(numel(params.x));
            obj.isCompositional = params.isCompositional;
            
        end
        
        % prior probability
        function [p, ln_p, grad_ln_p] = prob(obj,x)
            p = 1.0.*all(x == obj.mu);
            ln_p = log(p);
            grad_ln_p = zeros(numel(obj.mu),size(x,2));
        end
        
        % random sampler
        function x = rnd(obj,n)
            x = obj.mu .* ones(numel(obj.mu),n);        
        end
        
        % match first and second moments
        function obj = match(obj,x_mu,x_cov)
            obj.mu = x_mu;        
        end
        
    end

end