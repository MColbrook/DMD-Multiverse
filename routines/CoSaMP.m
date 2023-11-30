function [x,r,normR,residHist, errHist] = CoSaMP( A, b, k, errFcn, opts )
% x = CoSaMP( A, b, k )
%   uses the Compressive Sampling Matched Pursuit (CoSaMP) algorithm
%   (see Needell and Tropp's 2008 paper http://arxiv.org/abs/0803.2392 )
%   to estimate the solution to the equation
%       b = A*x     (or b = A*x + noise )
%   where there is prior information that x is sparse.
%
%   "A" may be a matrix, or it may be a cell array {Af,At}
%   where Af and At are function handles that compute the forward and transpose
%   multiplies, respectively. If the function handles are provided,
%   the the least-squares step is performed using LSQR (use could also use
%   CG on the normal equations, or other special variants).
%
% [x,r,normR,residHist,errHist] = CoSaMP( A, b, k, errFcn, opts )
%   is the full version.
% Outputs:
%   'x' is the k-sparse estimate of the unknown signal
%   'r' is the residual b - A*x
%   'normR' = norm(r)
%   'residHist'     is a vector with normR from every iteration
%   'errHist'       is a vector with the outout of errFcn from every iteration
%
% Inputs:
%   'A'     is the measurement matrix
%   'b'     is the vector of observations
%   'k'     is the estimate of the sparsity (you may wish to purposefully
%              over- or under-estimate the sparsity, depending on noise)
%              N.B. k < size(A,1) is necessary, otherwise we cannot
%                   solve the internal least-squares problem uniquely.
%   'errFcn'    (optional; set to [] to ignore) is a function handle
%              which will be used to calculate the error; the output
%              should be a scalar
%   'opts'  is a structure with more options, including:
%       .printEvery = is an integer which controls how often output is printed
%       .maxiter    = maximum number of iterations
%       .normTol    = desired size of norm(residual). This is also
%                       used to detect convergence when the residual
%                       has stopped decreasing in norm
%       .LSQR_tol   = when "A" is a set of function handles, this controls
%                       the tolerance in the iterative solver. For compatibility
%                       with older versions, the fieldname "cg_tol" is also OK.
%       .LSQR_maxit = maximum number of steps in the iterative solver. For compatibility
%                       with older versions, the fieldname "cg_maxit" is also OK.
%                       N.B. "CG" stands for conjugate gradient, but this code
%                            actually uses the LSQR solver.
%       .HSS        = if true, use the variant of CoSaMP that is similar to
%                       HHS Pursuit (see appendix A.2 of Needell/Tropp paper). Recommended.
%       .two_solves = if true, uses the variant of CoSaMP that re-solves
%                       on the support of size 'k' at every iteration
%                       (see appendix). This can be used with or without HSS variant.
%       .addK       = the number of new entries to add each time. By default
%                       this is 2*k (this was what is used in the paper).
%                       If you experience divergence, try reducing this.
%                       We recommend trying 1*k for most problems.
%       .support_tol = when adding the (addK) atoms with the largest
%                       correlations, the CoSaMP method specifies that you do
%                       not add the atoms if the correlation is exactly zero.
%                       In practice, it is better to not add the atoms of their
%                       correlation is nearly zero. "support_tol" controls
%                       what this 'nearly zero' number is, e.g. 1e-10.
%
%       Note that these field names are case sensitive!
%
%
% Stephen Becker, Aug 1 2011  srbecker@alumni.caltech.edu
%   Updated Dec 12 2012
%   See also OMP, test_OMP_and_CoSaMP

if nargin < 5, opts = []; end
if ~isempty(opts) && ~isstruct(opts)
    error('"opts" must be a structure');
end

function out = setOpts( field, default )
    if ~isfield( opts, field )
        opts.(field)    = default;
    end
    out = opts.(field);
end

printEvery  = setOpts( 'printEvery', 1000 );
maxiter     = setOpts( 'maxiter', 1000 );
normTol     = setOpts( 'normTol', 1e-10 );
cg_tol      = setOpts( 'cg_tol', 1e-6 );
cg_maxit    = setOpts( 'cg_maxit', 20 );
% Allow some synonyms
cg_tol      = setOpts( 'LSQR_tol', cg_tol );
cg_maxit    = setOpts( 'LSQR_maxit', cg_maxit );
HSS         = setOpts( 'HSS', false );
TWO_SOLVES  = setOpts( 'two_solves', false );
addK        = round(setOpts( 'addK', 2*k) );
support_tol = setOpts( 'support_tol', 1e-10 );


if nargin < 5 || isempty(printEvery)
    printEvery  = round(k,maxiter);
end

if nargin < 4
    errFcn = [];
elseif ~isempty(errFcn) && ~isa(errFcn,'function_handle')
    error('errFcn input must be a function handle (or leave the input empty)');
end

if iscell(A)
    LARGESCALE  = true;
    Af  = A{1};
    At  = A{2};     % we don't really need this...
else
    LARGESCALE  = false;
    Af  = @(x) A*x;
    At  = @(x) A'*x;
end


% -- Intitialize --
% start at x = 0, so r = b - A*x = b
r           = b;
Ar          = At(r);
N           = size(Ar,1);       % number of atoms
M           = size(r,1);        % size of atoms
if k > M/3
    error('K cannot be larger than the dimension of the atoms');
end
x           = zeros(N,1);
ind_k       = [];

% indx_set    = zeros(k,1);  % created on-the-fly
% A_T         = zeros(M,k);
residHist   = zeros(k,1);
errHist     = zeros(k,1);

fprintf('Iter,   |T|,  Resid');
if ~isempty(errFcn)
    fprintf(',   Error');
end
if LARGESCALE
    fprintf(',   LSQR iterations and norm(residual)');
end
fprintf('\n');

if LARGESCALE
      
%     if exist( 'lsqr_wrapper','file')
%         LSQR_ALG    = @lsqr_wrapper;
%         % This is Stephen's wrapper (designed to imitate Matlab's
%         %   syntax ) to the version of LSQR that you can get
%         % at http://www.stanford.edu/group/SOL/software/lsqr/matlab/lsqr.m
%         % This version of LSQR is better than Matlab's implementation.
%     elseif exist( 'lsqr','file' )
%         LSQR_ALG    = @lsqr; % Matlab's version
%     else
%         disp('You need to install LSQR! Download it from:');
%         disp('http://www.stanford.edu/group/SOL/software/lsqr/matlab/lsqr.m');
%         error('Need to have working copy of LSQR');
%     end

    % Unfortunately, the Stanford group's version doesn't make it easy
    %   to provide a starting value, so we will stick with Matlab's version.
    
    LSQR_ALG    = @(RHS,Afcn,x0) lsqr(Afcn,RHS,cg_tol,cg_maxit,[],[],x0 );  
end



for kk = 1:maxiter
    
    % -- Step 1: find new index and atom to add
    y_sort      = sort( abs(Ar),'descend');
    cutoff      = y_sort(addK); % addK is typically 2*k
    cutoff      = max( cutoff, support_tol );
    ind_new     = find( abs(Ar) >= cutoff );
    

    % -- Merge:
    T    = union( ind_new, ind_k );
    

    % -- Step 2: update residual
    if HSS
        RHS     = r; % where r = b - A*x, so we'll need to add in "x" later
        x_warmstart     = zeros(length(T),1);
    else
        RHS     = b;
        x_warmstart     = x(T);
    end
    
    
    % -- solve for x on the suppor set "T"
    if LARGESCALE
        % Could do CG on the normal equations...
        % Or, use LSQR:
%         x_T = LSQR_ALG(Afcn,RHS,cg_tol,cg_maxit);

        % use an initial guess to warm-start the solver
        Afcn = @(x,mode) partialA( N, T, Af, At, x, mode );
        [x_T,flag,relres,CGiter] = LSQR_ALG(RHS,Afcn,x_warmstart);  
        
    else
        x_T = A(:,T)\RHS;   % more efficient; equivalent to pinv when underdetermined.
%         x_T     = pinv( A(:,T) )*RHS;
    end
    
    
    if HSS
        % HSS variation of CoSaMP
        x_new       = zeros(N,1);
        x_new(T)    = x_T;
        x           = x + x_new;    % this is the key extra step in HSS
        cutoff      = findCutoff( x, k );
        x           = x.*( abs(x) >= cutoff );
        ind_k       = find(x);
        
        if TWO_SOLVES
            if LARGESCALE
                Afcn = @(x,mode) partialA( N, ind_k, Af, At, x, mode );
                [x_T2,flag,relres,CGiter] = LSQR_ALG(b,Afcn,x(ind_k)); % not using "r", just using "b"
            else
                x_T2     = A(:,ind_k)\b;
            end
            x( ind_k )  = x_T2;
        end
        
        % update r
        r_old   = r;
        r       = b - Af(x);
    else
        % Standard CoSaMP
        % Note: this is implemented *slightly* more efficiently
        %   that the HSS variation
    
        % Prune x to keep only "k" entries:
        cutoff  = findCutoff(x_T, k);
        Tk      = find( abs(x_T) >= cutoff );
        % This is assuming there are no ties. If there are,
        %    from a practical standpoint, it probably doesn't
        %    matter much what you do. So out of laziness, we don't worry about it.
        ind_k   = T(Tk);
        x       = 0*x;
        x( ind_k ) = x_T( Tk );
        
        
        if TWO_SOLVES
            if LARGESCALE
                Afcn = @(x,mode) partialA( N, ind_k, Af, At, x, mode );
                [x_T2,flag,relres,CGiter] = LSQR_ALG(b,Afcn,x(ind_k));
            else
                x_T2     = A(:,ind_k)\b;
            end
            x( ind_k )  = x_T2;
        end
        
        % Update x and r
        r_old   = r;
        if LARGESCALE
            r   = b - Af( x );
        else
            % don't do a full matrix-vector multiply, just use necessary columns
            r   = b - A(:,ind_k)*x_T( Tk );
        end
    end

    
    % -- Print some info --
    PRINT   = ( ~mod( kk, printEvery ) || kk == maxiter );
    normR   = norm(r);
    STOP    = false;
    if normR < normTol || norm( r - r_old ) < normTol
        STOP    = true;
        PRINT   = true;
    end
    
    
    if ~isempty(errFcn)
        er  = errFcn(x);
        errHist(kk)     = er;
    end
    if PRINT
        fprintf('%4d, %4d, %.2e', kk, length(T), normR);
        if ~isempty(errFcn)
            fprintf(', %.2e',er);
        end
        if LARGESCALE
            fprintf(',  %d,  %.2e', CGiter,relres);
        end
        fprintf('\n');
    end
        
    residHist(kk)   = normR;
    
    if STOP
        disp('Reached stopping criteria');
        break;
    end

    if kk < maxiter
        Ar  = At(r); % prepare for next round
    end
    
end

end % -- end of main function

function tau = findCutoff( x, k )
% finds the appropriate cutoff such that after hard-thresholding,
% "x" will be k-sparse
x   = sort( abs(x),'descend');
if k > length(x)
    tau = x(end)*.999;
else
    tau  = x(k);
end
end


% for LSQR, PCG, etc.
function z = partialA( N, T, Af, At, x, mode )
% Multiplies A_T(x) or A_T'(x) (the transpose)
%   where _T denotes restriction to the set "T"

switch lower(mode)
    case 'notransp'
        % zero-pad:
        y   = zeros(N,1);
        y(T)    = x;
        z       = Af(y);
        
    case 'transp'
        y   = At(x);
        % truncate:
        z   = y(T);
end
end % -- end of subfuction

