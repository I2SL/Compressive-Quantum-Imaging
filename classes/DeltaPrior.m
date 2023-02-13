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
            if numel(obj.mu)>1
                %p = double(ismember(x',obj.mu','rows')'); % exact vector equality
                p = double(vecnorm(x-obj.mu) < 1e-9); % tolerance vector equality
            else
                %p = double(x == obj.mu); % exact scalar equality
                p = double((x-obj.mu) < 1e-9);  % tolerance scalar equality
            end
            ln_p = log(p);
            grad_ln_p = zeros(numel(obj.mu),size(x,2));
        end
        
        % random sampler
        function x = rnd(obj,n)
            x = obj.mu .* ones(numel(obj.mu),n);        
        end
        
        % match first and second moments
        function obj = match(obj,x_mu,x_cov)
            %obj.mu = x_mu;  
            return
        end
        
    end

end