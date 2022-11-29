classdef MultilevelPrior < Prior
    properties
        % mean and covariance
        mu
        cov
        
        % compositional flag
        isCompositional
        
        % parameters for a multilevel distribution is just a collection of
        % prior objects
        P_objects
    end
    
    methods
        % constructor
        function obj = MultilevelPrior(params)
            % params should be a struct with fields:
            obj.P_objects = params.P_objects;               % a cell array of prior object instances
            
            % initialize isCompositional flag
            obj.isCompositional = 1;                        % a binary flag
            
            % construct mean, covariance and isCompositional parameters
            obj.mu = [];
            obj.cov = [];
            for k = 1:numel(obj.P_objects)
                P_k = obj.P_objects{k};
                obj.isCompositional = (obj.isCompositional & P_k.isCompositional); % multilevel Prior is compositional only if its constituents are all compositional
                obj.mu = vertcat(obj.mu, P_k.mu);
                obj.cov = blkdiag(obj.cov, P_k.cov);
            end
        end
        
        % prior probability
        function [p, ln_p, grad_ln_p] = prob(obj,x)
            p = 1;
            ln_p = 0;
            grad_ln_p = [];
            
            % initialize subvector start index
            k0 = 1;
            
            % apply the prob function for each constituent prior object
            for k = 1:numel(obj.P_objects)
                
                % get k'th prior object
                P_k = obj.P_objects{k};
                
                % number of parameters for k'th prior object
                n_k = numel(P_k.mu);
                
                % set subvector stop index for k'th prior object
                k1 = k0+(n_k - 1); 
                
                % call prob function for kth prior object and incorporate
                % results into multilevel outputs
                [p_k, ln_p_k, gln_p_k] = P_k.prob(x(k0:k1,:));
                p = p .* p_k;
                ln_p = ln_p + ln_p_k;
                grad_ln_p = vertcat(grad_ln_p, gln_p_k);
                
                % reset subvector start index
                k0 = k1 + 1; 
            end
            
        end
        
        % random sampler
        function x = rnd(obj,n)
            x = [];
            % apply the rnd function for each prior object to build up a
            % sample from the multilevel
            for k = 1:numel(obj.P_objects)
                % get k'th prior object
                P_k = obj.P_objects{k};
                
                % randomly sample k'th prior object and insert in random
                % sample
                x = vertcat(x, P_k.rnd(n));
            end
        end
        
        % match first and second moments
        function obj = match(obj, x_mu, x_cov)
            
                                
             % initialize subvector start index
            k0 = 1;
            
            % apply the prob function for each prior
            for k = 1:numel(obj.P_objects)
                
                % get k'th prior object
                P_k = obj.P_objects{k};
                
                % number of parameters for k'th prior object
                n_k = numel(P_k.mu);
                
                % set subvector stop index for k'th prior object
                k1 = k0+(n_k - 1); 
                
                % match moment for each constituent prior object
                obj.P_objects{k} = P_k.match(x_mu(k0:k1),x_cov(k0:k1,k0:k1));
                
                % reset subvector start index
                k0 = k1 + 1; 
            end

            
        end 
    end
    
end
