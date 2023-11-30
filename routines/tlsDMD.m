function [K,LAM,Phi] = tlsDMD(X,Y,r)
%%%%%%%%%%%%%%%%%%
% Applies total least-squares DMD
% References for algorithm:
% https://link.springer.com/article/10.1007/s00348-016-2127-7 and https://link.springer.com/article/10.1007/s00162-017-0432-2
% INPUTS: snapshot matrices X and Y, rank r
% OUTPUTS: Koopman matrix K, diagonal matrix of eigenvalues LAM, DMD modes
% Phi
%%%%%%%%%%%%%%%%%%

% First project onto POD modes
[U0,~,V] = svd([X,Y(:,end)],'econ');
r = min(min(rank(U0),floor(size(X,2)/2)),r);
U0 = U0(:,1:r);
Xtilde = U0'*X;
Ytilde = U0'*Y;

% Run the tls optimization
[U,~,~] = svd([Xtilde;Ytilde],'econ');

U11 = U(1:end/2,1:r);
U21 = U(end/2+1:end,1:r);

% Now extract the DMD ouputs
K = U21*pinv(U11); % DMD matrix
[W,LAM] = eig(K,'vector');

if (nargout > 2)
    Phi = Y*V*Sinv*W; % DMD modes
end
    
end


