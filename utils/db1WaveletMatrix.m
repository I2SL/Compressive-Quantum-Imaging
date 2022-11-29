function W = db1WaveletMatrix(max_depth,d)
% Creates the transform matrix where the columns are the vectorized 'db1'
% wavelets.

% assert max depth is a whole number
assert(max_depth == round(max_depth) && max_depth >= 0, 'decomposition depth is not a natural number');


dx = [ones(2^(max_depth-1),1);-ones(2^(max_depth-1),1)];

switch d
    case 1
        dim = [(2^max_depth),1]; % full-depth decomposition
        s = 1;                   % dc wavelet
        dy = 1;
    case 2
        dim = [1,1]*(2^max_depth); % full-depth decomposition
        s = ones(2^(max_depth),1); % dc wavelet
        dy = dx;
    otherwise
        error('signals with more than 2 dimensions are not supported')  
end

% assert requested decomposition level is valid 
assert(max_depth <= log2(max(dim)),'Wavelet decomposition level is deeper than allowed by image dimensions');


% construct db1 wavelet transform matrix
W = ones([dim,(2^max_depth)^d]); % automatically sets DC

k=2;
if max_depth > 0
    for level = (1:max_depth) - 1
        
        
        levelx = level;
        levely = 0;
        
        if d == 2
            levely = level;
        end
        
        % get scaling and new support dimensions
        scale = [2^-levelx,2^-levely];
        dimn = [dim(1)*scale(1),dim(2)*scale(2)];        
        
        for i = 0:(2^levelx -1)   
            for j = 0:(2^levely -1)
                
                % get shift
                shift = [i*dimn(1),j*dimn(2)];
                
                % get x diff vector
                dx_n = zeros(dim);
                dx_n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(dx*s',dimn,'nearest');
                W(:,:,k) =  dx_n; k = k+1;


                % get y diff and xy diff vectors
                if d == 2
                    dy_n = zeros(dim);
                    dxdy_n = zeros(dim);
                    dy_n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(s*dy',dimn,'nearest');
                    dxdy_n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(dx*dy',dimn,'nearest');
                    W(:,:,k) = dy_n; k = k+1;
                    W(:,:,k) = dxdy_n; k = k+1;
                end
            end
        end
    end
end



% vectorize all of the 2D wavelets
W = squeeze(reshape(W,[prod(dim),1,(2^max_depth)^d]));

% make the matrix orthogonal
W = W ./ sqrt((sum(abs(W).^2,1)));
end


%{
% append the DC wavelet 
W(:,:,1) = s;

% The difference elemental db1 wavelets (diffx, diffy, diffxy)
w1 = vertcat(ones(dim(1)/2,dim(2)),-ones(dim(1)/2,dim(2)));                 % [1,  1 ; -1, -1]
w2 = horzcat(ones(dim(1),dim(2)/2),-ones(dim(1),dim(2)/2));                 % [-1, 1 ; -1,  1]
w3 = w1.*w2;                                                                % [-1, 1 ;  1, -1]

k = 1;

% append the higher-order difference wavelets
if max_depth > 0      
    for level = (1:max_depth)-1
        for i = 0:(2^level -1)
            for j = 0:(2^level -1)
                w1n = zeros(dim);
                w2n = zeros(dim);
                w3n = zeros(dim);

                scale = 2^-level;
                dimn = dim * scale;
                shift = [i*dimn(1),j*dimn(2)];

                w1n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(w1,scale,'nearest');
                w2n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(w2,scale,'nearest');
                w3n(shift(1)+(1:dimn(1)),shift(2)+(1:dimn(2))) = imresize(w3,scale,'nearest');

                W(:,:,k+1) = w1n;
                W(:,:,k+2) = w2n;
                W(:,:,k+3) = w3n;

                k = k+3;
            end
        end
    end
end


% vectorize all of the 2D wavelets
W = squeeze(reshape(W,[prod(dim),1,(2^max_depth)^d]));
end
%}
