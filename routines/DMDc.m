function [A,B,LAM,Phi] = DMDc(X,Y,Upsilon,p,r)
%%%%%%%%%%%%%%%%%%
% Applies DMD with control
% INPUTS: snapshot matrices X, Y and Upsilon, rank r
% OUTPUTS: Matrices A and B, diagonal matrix of eigenvalues LAM, DMD modes
% Phi
%%%%%%%%%%%%%%%%%%

[U,S,V] = svd([X;Upsilon],'econ');
p = min(rank(S),p);
U = U(:,1:p); V = V(:,1:p); S = S(1:p,1:p); Sinv = diag(1./diag(S));

U1 = U(1:size(X,1),:);
U2 = U(size(X,1)+1:end,:);

[Uhat,~,~] = svd(Y,'econ');
r = min(rank(Uhat),r);
Uhat = Uhat(:,1:r); 

A = (Uhat')*Y*V*Sinv*(U1')*Uhat;
B = (Uhat')*Y*V*Sinv*(U2');

[W,LAM] = eig(A,'vector');

if (nargout > 3)
    Phi = Y*V*Sinv*(U1')*Uhat*W; % DMD modes
end
    
end