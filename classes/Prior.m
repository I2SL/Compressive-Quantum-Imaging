classdef (Abstract) Prior
    properties (Abstract)
        % mean and covariance
        mu
        cov
        
        % A flag indicating if the prior is compositional.
        % These are priors that automatically satisfy the positivity constraint.
        % https://en.wikipedia.org/wiki/Compositional_data
        isCompositional
    end
    methods (Abstract)
        % prior probability
        prob(obj,x)
        % random sampler
        rnd(obj,n)
        % match first and second moments
        match(obj,x_mu,x_cov)
    end
end