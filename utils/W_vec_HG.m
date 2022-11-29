function  W_vec = W_vec_HG(W,rl,n_modes,n_samples,d)
% Computes the wavelet operators in the HG representation for a given
% wavelet transform. These operators are related to a wavelet transform of
% the system PSF, which is assumed to be a non-skew 2D gaussian with
% width sigma_x = sigma_y.
% ----------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------
% W             : wavelet transform matrix. Each column in W is a wavelet
%                 vector from an orthogonal and complete wavelet basis.
% rl            : the rayleigh limit (as a fraction of the signal support)
% n_modes       : Maximum order of 1D Hermite-Gauss modes used for the matrix
%                 representation. The total number of 2D modes used is thus n(n+1)/2
% n_samples     : number of samples per dimension to use for each HG mode
% d             : dimension of signal d = 1,2

% ----------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------
% W_vec         : a matrix stack of wavelet operators


% upsample the wavelet transform matrix
N_samples = n_samples^d;
N_params = size(W,2);
W = imresize(W,[N_samples,N_params],'nearest')*size(W,1)/N_samples;

% Gaussian PSF width
sigma = rl;

% discretize the image plane
x = linspace(-0.5,+0.5,n_samples);
dx = x(2)-x(1);                     % differential element
xx = x/2/sigma;                     % dimensionless coordinates

% get the total number of modes
switch d
    case 1
        N_modes = n_modes;
        yy = 0;        
    case 2
        N_modes = n_modes*(n_modes+1)/2;
        yy = xx;
    otherwise
      error('signals with more than 2 dimensions are not supported')  
end
[XX,YY] = meshgrid(xx,yy);
        
count = 1;
% build up index
for i = 1:n_modes
    for j = 1:(i^(d-1))

        HG_proj(count).ind_x = i-j;
        HG_proj(count).ind_y = i-HG_proj(count).ind_x-1;  

        count = count + 1;

    end
end

W_vec = zeros(N_modes,N_modes,N_params);

for i = 1:N_modes
    for j = 1:N_modes
        p = HG_proj(i).ind_x;
        q = HG_proj(i).ind_y;
        m = HG_proj(j).ind_x;
        n = HG_proj(j).ind_y;
        
        % g(X,Y) = <HG_pq|Psi(X,Y)><Psi(X,Y)|HG_mn>
        g = XX.^(p+m) .* YY.^(q+n) .* ...
            exp(- (XX.^2 + YY.^2)) * ...
            1/sqrt(factorial(p)*factorial(q) ...
            *factorial(m)*factorial(n));
        
       % vectorize g
       g = reshape(g,[numel(xx)*numel(yy),1]);
       
       % wavelet transform g(X,Y)       
       w = W'*g; % equivalent to W\g if W is an orthonormal matrix
       
       % assign wavelet coefficients from each wavelet to matrix element in
       % stack
       W_vec(i,j,:) = w;        
    end
end

end



%{
% -------------------------------------------------------------------
% Gaussian PSF width
sigma = rl;
sigma_x = sigma;    
sigma_y = sigma;

% Discretize image plane coordinates system
x = linspace(-0.5,+0.5,n_samples);
y = x;
[X,Y] = meshgrid(x,y);

% dimensionless coordinates
XX = (X/2/sigma_x); 
YY = (Y/2/sigma_y);

% upsample the wavelet transform matrix
N_samples = n_samples^2;
N_params = size(W,2);
W = imresize(W,[N_samples,N_params],'nearest');

% get number of modes
N_modes = n_modes*(n_modes+1)/2;
A_vec = zeros(N_modes,N_modes,N_params); 

count = 1;
% build up index
for i = 1:n_modes
    for j = 1:i

        HG_proj(count).ind_x = i-j;
        HG_proj(count).ind_y = i-HG_proj(count).ind_x-1;  
        
        count = count + 1;

    end
end




for i = 1:N_modes
    for j = 1:N_modes
        p = HG_proj(i).ind_x;
        q = HG_proj(i).ind_y;
        m = HG_proj(j).ind_x;
        n = HG_proj(j).ind_y;
        
        % g(X,Y) = <HG_pq|Psi(X,Y)><Psi(X,Y)|HG_mn>
        g = XX.^(p+m) .* YY.^(q+n) .* ...
            exp(- (XX.^2 + YY.^2)) * ...
            1/sqrt(factorial(p)*factorial(q) ...
            *factorial(m)*factorial(n));
        
       % vectorize g
       g = reshape(g,[numel(x)*numel(y),1]);
       
       % wavelet transform g(X,Y)       
       w = W'*g; % equivalent to W\g if W is an orthonormal matrix
       
       % assign wavelet coefficients from each wavelet to matrix element in
       % stack
       A_vec(i,j,:) = w;        
    end
end
%}
