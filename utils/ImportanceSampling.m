function [x_mu, x_cov, valid_count] = ImportanceSampling(post_fn, ref_fn, ref_rnd, N_samples, use_log_shift)
    % post_pdf -  a function handle to the propritional pdf of the
    %             posterior
    % ref_rnd   - a random sample generator from the reference distribution.
    %             Outputs an MxN matrix where M is the dimensionality of
    %             the samples and N is the number of samples
    % ref_pdf   - a function handle to the pdf of reference distribution
    % N_samples - the number of samples used to compute the posterior
    %             moments
    %
    % References: 
    % https://dept.stat.lsa.umich.edu/~jasoneg/Stat406/lab7.pdf

    % draw samples from the reference distribution
    samples = ref_rnd(N_samples);
    
    % evaluate the reference function at the samples
    [ref_pdf,log_ref_pdf,~] = ref_fn(samples);
    
    % evaluate the posterior function at the drawn samples
    [post_pdf, log_post_pdf,~] = post_fn(samples);
    
    % apply log shift to relative probabilities to increase the number of
    % valid samples
    if use_log_shift
        % compute log shift amount that ensures pdf_ratio =/= inf
        s_i  = log(realmax/N_samples) + log_ref_pdf - log_post_pdf;
        s = min(s_i);
        
        % apply shift
        log_post_pdf = log_post_pdf + s;
        
        % compute the ratio of the pdfs for each sample with log shift
        pdf_ratio = exp(log_post_pdf - log_ref_pdf);
    else
        % compute the ratio of pdfs for each sample
        pdf_ratio = post_pdf ./ ref_pdf;
    end
    
    % remove invalid samples
    pdf_ratio(isnan(pdf_ratio)) = 0;
    pdf_ratio(isinf(pdf_ratio)) = 0; % maybe set this to a large number
    
    % get index of the non-zero pdf ratio elements
    valid_idx = pdf_ratio>0;                % valid samples binary index
    pdf_ratio = pdf_ratio(valid_idx);       % reject pdf ratio elements that are invalid
    samples = samples(:,valid_idx);         % reject samples that are invalid
        
    % count the number of valid samples  (samples with numerically stable probability) 
    % after the log shift
    valid_count = nnz(valid_idx);
    disp(['Valid Samples: ', num2str(valid_count)])
    
    assert(valid_count > 0,'Importance sampling did not converge.')
    
    % normalizing constant
    C = sum(pdf_ratio);
    
    % measure for good reference distribution
    ESS = sqrt(mean((pdf_ratio./mean(pdf_ratio) - 1).^2)); % should be less than 5 for a good estimator
    
    % estimate the first moment of the posterior (expected value)
    x_mu = sum(pdf_ratio .* samples,2) / C;
    
    % estimate the second moment of the posterior
    sv = reshape(samples,[size(samples,1),1,size(samples,2)]);
    sh = reshape(samples,[1,size(samples,1),size(samples,2)]);
    pr = reshape(pdf_ratio,[1,1,numel(pdf_ratio)]);
    xx = sum(pr .* sv .* sh ,3) / C;
       
    % estimate the covariance of the posterior
    x_cov = xx - x_mu*x_mu';
    
    
    % make covariance positive semidefinite
    epsilon = 1e-10;
    ev = eig(x_cov);
    if ~all(ev >0)
        x_cov = x_cov + (epsilon - min(ev))*eye(size(x_cov));
    end
    
    
    
    %{
    % approximate the covariance of the posterior
    delta = samples - x_mu;
    
    % (compatible with Matlab 2017)
    n = numel(x_mu);
    x_cov = zeros([n,n]);
    for i = 1:valid_count
        x_cov = x_cov + pdf_ratio(i) * (delta(:,i) * delta(:,i)');
    end
    
    x_cov = x_cov;
    %}
    
    %{
    % visualize the log posterior density in 3D
    show_samples = post_pdf(valid_idx);
    figure
    av = [1,2,3];   % indices of a parameters we'd like to vizualize
    sz = 5*(show_samples-min(show_samples))/(max(show_samples)-min(show_samples));    % normalized size
    scatter3(samples(av(1),:),samples(av(2),:),samples(av(3),:),sz+1e-6,'filled');
    hold on
    scatter3(x_mu(av(1)),x_mu(av(2)),x_mu(av(3)),'filled','r')
    hold off
    %}
    
end