classdef Posterior
   properties
       num_samples = 1e5;
       K = 1;       
       M_set = [1];
       prior_obj
       ref_obj
       log_shift = 1;
       method = 'cumulative' %{'cumulative','conjugate','K_informative'}
       
       
   end
   methods
       function obj = Posterior(prior_obj, ref_obj)
           % set object properties
           obj.prior_obj = prior_obj;
           obj.ref_obj = ref_obj;
       end
       
       function [obj,x_mu,x_cov] = PosteriorMoments(obj,l_seq,B_seq,W_vec,mu_seq,cov_seq)
          
          % get the iteration
          iter = size(l_seq,2);
          
          if iter > 1
              % set the subset of measurement indices based on the method
              switch obj.method
                  case 'cumulative'
                      obj.ref_obj = obj.ref_obj.match(mu_seq(:,iter),cov_seq(:,:,iter));
                      obj.M_set = 1:iter;

                  case 'conjugate'
                      obj.ref_obj = obj.ref_obj.match(mu_seq(:,iter),cov_seq(:,:,iter));
                      obj.prior_obj = obj.prior_obj.match(mu_seq(:,iter),cov_seq(:,:,iter));
                      obj.M_set = iter;

                  case 'K_informative'  
                      %obj.ref_obj = obj.ref_obj.match(mu_seq(:,iter),cov_seq(:,:,iter));
                      %obj.prior_obj = obj.prior_obj.match(mu_seq(:,iter),cov_seq(:,:,iter));
                      
                      
                      % greedily determine the best subset of measurements
                      % to use when calculating the cumulative likelihood
                      if iter <= obj.K
                          obj.M_set = 1:iter;
                      else 
                          k_set = ~eye(obj.K);
                          Q_Mk = zeros(obj.K,1);
                            
                          % determine the merit of all possible subsets
                          % involving the new measurment
                          for k = 1:obj.K
                              M_set_k = [obj.M_set(k_set(k,:)),iter]; 
                              Q_Mk(k) = obj.Merit(l_seq(:,M_set_k), B_seq(:,:,M_set_k), W_vec);                         
                          end
                          [~,k_opt] = max(Q_Mk);
                          obj.M_set = [obj.M_set(k_set(k_opt,:)), iter];
                      end
              end
          end
          
          B_set = B_seq(:,:,obj.M_set);
          l_set = l_seq(:,obj.M_set);
          
          % function handles for Importance sampling
          ref_fn = @(x) obj.ref_obj.prob(x);
          ref_rnd = @(n) obj.ref_obj.rnd(n);
          post_fn = @(x) obj.prob(x, l_set, B_set, W_vec);
          
          [x_mu, x_cov] = ImportanceSampling(post_fn, ref_fn, ref_rnd, obj.num_samples, obj.log_shift);
                  
       end
       
       function [p, ln_p, grad_ln_p] = prob(obj,x,l_set,B_set,W_vec)
           % get likelihood and prior
           [like_p, like_ln, like_grad_ln] = likelihood(x,l_set,B_set,W_vec);
           [prior_p, prior_ln, prior_grad_ln] = obj.prior_obj.prob(x);
           
           % combine outputs for posterior
           p = like_p .* prior_p;
           ln_p = like_ln + prior_ln;
           grad_ln_p = like_grad_ln + prior_grad_ln;           
       end
       
       function merit = Merit(obj, l_set, B_set, W_vec)
           
          % function handles for Importance sampling
          ref_fn = @(x) obj.ref_obj.prob(x);
          ref_rnd = @(n) obj.ref_obj.rnd(n);
          post_fn = @(x) obj.prob(x, l_set, B_set, W_vec);
          
          [x_mu, x_cov] = ImportanceSampling(post_fn, ref_fn, ref_rnd, obj.num_samples, obj.log_shift);
          
          % compute merit of the measurement set 
          [~, E_Q] = JointParameterEstimator(W_vec, x_mu, x_cov);
          merit = - min( eig(E_Q) );
       end
      
      
   end
end