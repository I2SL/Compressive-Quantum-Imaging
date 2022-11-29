function V = V_matrix(W,M)
    % get polytope constraint matrix (constraint 2: non-negativity)
    HH = W*M;
    H = -HH(:,2:end);                     % half-space matrix
    c = HH(:,1);                          % half-space vector 
    V = lcon2vert(H,c)';                  % convert H-representation (H*a <= c) to V-representation
    V = padarray(V,[1,0],1,'pre');        % augment polytope vertex matrix: (aa = V*x where x is from the N-params Simplex)
end