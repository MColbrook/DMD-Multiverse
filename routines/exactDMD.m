function [K,LAM,Phi] = exactDMD(X,Y,r)
%%%%%%%%%%%%%%%%%%
% Applies exact DMD
% INPUTS: snapshot matrices X and Y, rank r
% OUTPUTS: Koopman matrix K, diagonal matrix of eigenvalues LAM, DMD modes
% Phi
%%%%%%%%%%%%%%%%%%

[U,S,V] = svd(X,'econ');
r = min(rank(S),r);
U = U(:,1:r); V = V(:,1:r); S = S(1:r,1:r); Sinv = diag(1./diag(S));
K = (U')*Y*V*Sinv; % DMD matrix
[W,LAM] = eig(K,'vector');

if (nargout > 2)
    Phi = Y*V*Sinv*W; % DMD modes
end
    
end