classdef SimulatedMeasurement < Measurement
    properties
        rho         % ground truth state density
        n_pho       % number of photons collected per call to measurement
        m_available % number of calls to measurement available
    end
    
    methods
        % constructor
        function obj = SimulatedMeasurement(rho,n_pho,max_iter)
            % checks
            assert(abs(trace(rho) - 1) < 1e-3, 'the density operator is not trace 1');
            assert(rem(n_pho,1)==0 & n_pho > 0, 'number of photons is not a positive integer');
            assert(rem(max_iter,1)==0 & max_iter > 0, 'max number of available measurements is not a positive integer');
            
            % set class properties
            obj.stop_flag = max_iter == 0;
            obj.rho = rho;
            obj.n_pho = n_pho;
            obj.m_available = max_iter;   

        end
        
        % measurement
        function [obj,l_vec] = MakeMeasurement(obj,B_gamma)
        % Simulates the result of a photon counting measurement with the
        % measurement operator B_gamma (the joint parameter estimator).
        % The interrogated mixed state is the ground-truth density operator obj.rho.
        % Each element of the output is the number of photons observed in 
        % each eigenstate of B_gamma. 
        %
        % --------
        % Inputs:
        % --------
        % B_gamma - a Hermitian matrix corresponding to an observable. In this case
        %           it is the joint parameter estimator.
        % --------
        % Outputs:
        % --------
        % l_vec - a column vector containing the number of photons observed in different 
        %          eigenstates of B_gamma.
        %%
            
            % return no photons if the stop flag is set
            if obj.stop_flag
                l_vec = [];
            else
                % Get measurement projectors (eigenstates) of the
                % the joint parameter operator
                [V,~] = eig(B_gamma);

                % From the measurement projectors, get the probability distribution over
                % eigenvalues of B_gamma of observing B_gamma. This is the true likelihood
                % p(l_vec | rho)

                % Two requirements for the probability to be valid
                %-----assert(all(V'*V == eye(size(V,1))));
                %-----assert(abs(trace(gt_rho)-1)<1e-10,'Density operator is not trace 1.')

                % container for photon detection probabilities in each of the eigenstates
                % of the measurement operator
                D = size(B_gamma,1); % dimensionality of Hilbert space
                p_eig = zeros(D+1,1);

                % get detection probabilities
                p_eig(1:D) = sum(V'*obj.rho.*V',2); % equivalent to diag(V'*gt_rho*V) but faster since off-diagonal elements are not computed

                % throw a warning if some of the probabilities are negative (numerical
                % instabilitiy
                neg_prob = p_eig<0;
                if any(neg_prob)
                    warning('Negative probabilties computed for modal measurement outcomes in measurement simulation. Setting negative probabilities to zero. Consider up-sampling the image plane.')
                    p_eig(neg_prob) = 0;
                end

                % compute the detection probability for modes beyond the N_HG_modes
                residual = 1-sum(p_eig,1);
                residual(residual < 0) = 0; % no residual if probabilities already sum to 1 (or more)

                % include the residual probability 
                p_eig(D+1) = residual;

                % normalize probability
                p_eig = p_eig/sum(p_eig);

                %% Generate samples from the likelihood
                % sample the modal indices in which each detected photon is found from the discrete distribution
                % p_eig associated with the probability of single-photon detections in the measurement eigenstates
                measurement = datasample(1:numel(p_eig),obj.n_pho,'Weights',p_eig)';
    
                % l_vec contains the total photon count in each mode
                [gc,gs] = groupcounts(measurement);
                l_vec = zeros(1,numel(p_eig));
                l_vec(gs) = gc;
                
                
                %{
                ind = datasample(1:numel(p_eig),obj.n_pho,'weights',p_eig)';

                % l_vec contains the count of the number of times each outcome appeared
                l_vec = accumarray(ind, 1);

                % fill in the cases where no counts appeared
                min_ind = min(ind);
                max_ind = max(ind);
                l_vec = [zeros(1-min_ind,1);l_vec;zeros(numel(p_eig)-max_ind,1)];
                %}
                
                % decrement the number of measurement available and trigger
                % stop flag if no measurements are left.
                obj.m_available = obj.m_available -1;
                if obj.m_available == 0
                    obj.stop_flag = 1;
                end
            end
            
        end
    end
end