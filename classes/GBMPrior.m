classdef GBMPrior < Prior
    properties
        
        % mean and covariance
        mu
        cov
        
        % compositional flag
        isCompositional = 0
 
        % parameters for a Gaussian Binomial Mixture distribution
        q
        mu0
        mu1
        cov0
        cov1
    end
    
    methods
        % constructor
        function obj = GBMPrior(params)
            % params should be a struct with fields:
            obj.q = params.q;
            obj.mu0 = params.mu0;
            obj.mu1 = params.mu1;
            obj.cov0 = params.cov0;
            obj.cov1 = params.cov1;
            
            % assign mean and covariance
            obj.mu = (1-obj.q)*obj.mu0 + obj.q*obj.mu1;
            obj.cov = (1-obj.q) * (obj.cov0 + obj.mu0*obj.mu0') + ...
                         obj.q  * (obj.cov1 + obj.mu1*obj.mu1') - ...
                         obj.mu * obj.mu';
        end
        
        % prior probability
        function [p, ln_p, grad_ln_p] = prob(obj,x)
            % prior
            p = ((1-obj.q)*mvnpdf(x',obj.mu0',obj.cov0) +...
                obj.q*mvnpdf(x',obj.mu1',obj.cov1))';

            % log prior
            ln_p = log(p);

            % gradient of prior
            grad_p = -(1-obj.q) * mvnpdf(x',obj.mu0',obj.cov0)'.*(obj.cov0\(x-obj.mu0))...
                     - obj.q * mvnpdf(x',obj.mu1',obj.cov1)'.*(obj.cov1\(x-obj.mu1)); % covariance matrix inverses may cause instability 

            % gradient of log prior
            grad_ln_p = grad_p ./ p; 
        end
        
        % random sampler
        function x = rnd(obj,n)
            
            % produces n random samples GBM vectors
            coinflips = binornd(1,obj.q,[n,length(obj.mu0)]);
            mvn0 = mvnrnd(obj.mu0',obj.cov0,n); 
            mvn1 = mvnrnd(obj.mu1',obj.cov1,n);
            x = (~coinflips.*mvn0 + coinflips.*mvn1)';
            
        end
        
        % match first and second moments
        function obj = match(obj,x_mu,x_cov) 
            % match moments with non-sparse normal distribution.
            obj.cov1 = x_cov;
            obj.mu1 = x_mu;
            %{
            % only updates parameters for N1 (not N0) to match moments
            obj.cov1 = ( (x_cov + x_mu*x_mu') - (1-obj.q)*(obj.cov0 + obj.mu0*obj.mu0')) / (obj.q) - obj.mu1 * obj.mu1';       
            obj.mu1 = (x_mu - (1-obj.q)*obj.mu0)/obj.q;
            %}
        end
    end
        
    
end
