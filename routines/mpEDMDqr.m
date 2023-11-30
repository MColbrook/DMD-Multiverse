function [mpK,mpV,mpD] = mpEDMDqr(PX,PY,W)
% mpEDMD algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   PX: dictionary evaluated at snapshots,
%           PY: dictionary evaluated at snapshots one time step later
%           W:  vector of weights for the quadrature
% OUTPUTS:   Koopman matrix mpK and its eigendecomposition [mpV,mpD]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLEASE CITE:  Colbrook, Matthew J. "The mpEDMD algorithm for data-driven computations of measure-preserving dynamical systems." 
%                                    SIAM Journal on Numerical Analysis 61.3 (2023): 1585-1608.
% Author and copyright:  Matthew Colbrook, https://www.damtp.cam.ac.uk/user/mjc249/home.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = W(:);
[Q,R] = qr(sqrt(W).*PX,"econ") ;
T = (R')\(PY'*(sqrt(W).*Q));
[U,~,V] = svd(T);
[mpV,mpD] = schur(V*U','complex'); % Schur decomposition used to ensure orthonormal basis

mpV = R\mpV;
mpK = (R\(V*U'))*R;
mpD = diag(diag(mpD));

end




