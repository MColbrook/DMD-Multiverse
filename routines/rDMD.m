function [K,LAM,tm,Phi] = rDMD(X,Y,r,p)
Z = [X,Y(:,end)];

tic
Q = randQB(Z,r,p);

X = Q'*Z;
Y = X(:,2:end);
X = X(:,1:end-1);

[U,S,V] = svd(X,'econ');
r = min(rank(S),r);
C = Y*V(:,1:r)*diag(1./diag(S(1:r,1:r)));
K = (U(:,1:r))'*C; % DMD matrix
[W,LAM] = eig(K,'vector');

if (nargout > 3)
    Phi = Q*C*W; % DMD modes
end
tm = toc;


end



function Q = randQB(X,k,p)

OM = randn(size(X,2),k+p);
Y = X*OM;

% for j = 1:q
%     [Q,~] = qr(Y,"econ");
%     [Z,~] = qr(X'*Q,"econ");
%     Y = X*Z;
% end

[Q,~]=qr(Y,"econ");

end

