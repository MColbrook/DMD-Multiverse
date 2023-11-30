function [K,LAM,Phi] = fbDMD(X,Y,r)
%%%%%%%%%%%%%%%%%%
% Applies forward-backward DMD
% Reference for algorithm: https://link.springer.com/article/10.1007/s00348-016-2127-7
% INPUTS: snapshot matrices X and Y, rank r
% OUTPUTS: Koopman matrix K, diagonal matrix of eigenvalues LAM, DMD modes
% Phi
%%%%%%%%%%%%%%%%%%

[U,S,V] = svd([X,Y(:,end)],'econ');
r = min(rank(U),r);
U = U(:,1:r); V = V(:,1:r); S = S(1:r,1:r); Sinv = diag(1./diag(S));
X2 = U'*X; Y2 = U'*Y;

[U1,S1,V1] = svd(X2,'econ');
[U2,S2,V2] = svd(Y2,'econ');
r = min(min(rank(S1),rank(S2)),r);

U1 = U1(:,1:r); V1 = V1(:,1:r); S1 = S1(1:r,1:r); Sinv1 = diag(1./diag(S1));
K1 = (U1')*Y2*V1*Sinv1; % forward DMD matrix

U2 = U2(:,1:r); V2 = V2(:,1:r); S2 = S2(1:r,1:r); Sinv2 = diag(1./diag(S2));
K2 = (U2')*X2*V2*Sinv2; % backward DMD matrix

[wf,ef] = eig(K1);
wf = Y2*(V1*(S1\wf));

[wb,eb] = eig(K2);
wb = X2*(V2*(S2\wb));

K1 = wf*ef*pinv(wf);
K2 = wb*eb*pinv(wb);

% atilde = sqrtm(fatilde*pinv(batilde));
% [atilde,errtemp] = besterrnbyn(fatilde,atilde,ifoverride);
% 
% [w,e] = eig(atilde);
% e = diag(e);


K = sqrtm(K1*pinv(K2));
% [K,~] = besterrnbyn(K1,K);


[W,LAM] = eig(K,'vector');

% choose the signs of the eigenvalues
Kfit = diag(W\(K1*W));
for jj = 1:size(Kfit)
    if abs(Kfit(jj)+LAM(jj))<abs(Kfit(jj)-LAM(jj))
        LAM(jj)=-LAM(jj);
    end
end

[~,I]=sort(abs(1-LAM),'ascend');
W =W(:,I); LAM =LAM(I);

K = W*diag(LAM)/W;

if (nargout > 2)
    Phi = Y*V*Sinv*W; % DMD modes
end
    
end



function [Cbest,besterr] = besterrnbyn(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find modification of B which is closest to
% A. We assume that B is computed as a square root.
% The nonuniqueness of matrix square roots
% is up to +/- times each eigenvalue (2^n 
% possibilities for an n x n matrix)
%
% Warning: this is a slow routine! It is exponential
% in n and should not be called for large n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n] = size(A);

if n > 25
    error('this will take forever. bomb')
end

C = B;
besterr = norm(A-C,'fro');
Cbest = C;

[w,e] = eig(B);
de = diag(e);
pw = pinv(w);

for i = 1:2^n-1
    v = (-1).^str2num(dec2bin(i,n)');
    C = w*diag(de.*v)*pw;
    err = norm(A-C,'fro');
    if (err < besterr)
        besterr = err;
        Cbest = C;
    end
end

end


