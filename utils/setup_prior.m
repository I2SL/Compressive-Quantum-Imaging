function prior_obj = setup_prior(prior_name)
% ammend this switch-case list to incorporate desired priors
switch prior_name
    
    case 'dirichlet'
        params.alpha = [];
        
        prior_obj = DirichletPrior(params);

    case 'mvn'
        params.mu = [];
        params.cov = [];
        
        prior_obj = MVNPrior(params);
        
    case 'gbm'
         params.q = [];
         params.cov0 = [];
         params.cov1 = [];
         params.mu0 = [];
         params.mu1 = [];

         prior_obj = GBMPrior(params);
        
    case 'cifar10gbm'
        
        % load cifar10 EM parameters
        load GBM_params_level.mat q_level cov0_level cov1_level

        % instantiate a multilevel GBM prior
        depth = 2;
        params.P_objects = cell(1+depth,1);

        % level 0
        prms_k = struct();
        prms_k.isCompositional = 0;
        prms_k.x = 1;
        params.P_objects{1} = DeltaPrior(prms_k);

        % level > 0
        for k = 1:depth

            % GBM parameter struct
            prms_k = struct();
            prms_k.q = q_level(k);
            prms_k.cov0 = cov0_level{k};
            prms_k.cov1 = cov1_level{k};
            prms_k.mu0 = zeros(size(prms_k.cov0,1),1);
            prms_k.mu1 = prms_k.mu0;

            % add GBM object to object collection
            params.P_objects{k+1} = GBMPrior(prms_k);
        end


        prior_obj = MultilevelPrior(params);
end

end



