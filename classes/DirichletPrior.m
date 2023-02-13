classdef DirichletPrior < Prior
    properties
        
        % mean and covariance
        mu
        cov
        
        % compositional flag
        isCompositional = 1
        
        % parameters for the Dirichlet distribution
        alpha
    end
    
    methods
         % constructor
        function obj = DirichletPrior(params)
            % checks on construction
            assert(all(params.alpha > 0),'dirichlet hyperparameters are not postive')
            
            % params should be a struct with fields:
            obj.alpha = params.alpha;
            
            % assign mean and covariance
            a0 = sum(obj.alpha);
            obj.mu = obj.alpha/a0;
            obj.cov = (diag(obj.alpha)/a0 - (obj.alpha * obj.alpha')/(a0^2))/(a0+1);
        end
        
        % prior probability
        function [p, ln_p, grad_ln_p] = prob(obj,x)
            
            % check input
            assert(all(x(:)>=0) && all(sum(x,1)-1<1e-5),'dirichlet variable lies outside simplex');
    
            % prior
            p = prod(gamma(obj.alpha)) / gamma(sum(obj.alpha)) .* prod(x.^(obj.alpha-1),1);

            % log prior
            ln_p = sum((obj.alpha-1).*log(x),1) + sum(log(gamma(obj.alpha))) - log(gamma(sum(obj.alpha)));

            % gradient of log prior
            grad_ln_p = obj.alpha ./ x;
            
        end
        
        % random sampler
        function x = rnd(obj,n)
            % returns n samples of the dirichlet distribution parametrized by alpha
            
            % Method for sampling from dirichlet using gamma distributions:
            % https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
            y = gamrnd(repmat(obj.alpha,[1,n]),1);
            x = y./sum(y,1);
        end
        
        % match first and second moments
        function obj = match(obj,x_mu,x_cov)
            d_cov = diag(x_mu) - x_mu*x_mu'; % dirichlet distribution covariance structure
            % determine optimal scaling parameter
            a0 = (x_cov(:)'*d_cov(:)/(x_cov(:)'*x_cov(:))) - 1;
            % reset the hyperparameters for the dirichlet distribution
            obj.alpha = a0*x_mu;
        end
    end
end