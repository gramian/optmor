function XP = optmor(f,g,s,t,r,q,co,nf,ut,x0,yd)
% optmor (Version 2.1)
% by Christian Himpe, 2013-2015 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNTAX:
%    W = optmor(f,g,s,t,r,q,[co],[nf],[ut],[x0],[yd]);
%
% SUMMARY:
%    optmor - optimization-based model order reduction,
%    computation of empirical gramians for model reduction,
%    system identification and uncertainty quantification.
%    Enables gramian-based nonlinear model order reduction.
%    Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%   (func handle)  f - system function handle; signature: xdot = f(x,u,p)
%   (func handle)  g - output function handle; signature:    y = g(x,u,p)
%        (vector)  s - system dimensions [inputs,states,outputs]
%        (vector)  t - time discretization [start,step,stop]
%        (scalar)  r - reduced order or error threshold
%        (vector)  q - parameter
%        (matrix,vector,scalar)  [co = 1] - covariance matrix
%        (vector,scalar) [nf = 0] - options, 10 components:
%            + Optimization Algorithm: fminunc(0), fminsearch(1), Custom Optimizer(-1)
%            + Lasso Regularization Weight: default(0)
%            + Tikhonov Regularization Weight: default(0)
%            + Data-Driven Regularization Weight: default(0)
%            + Monte-Carlo Basis Enrichment (Base Size)
%            + Random Seed: none(0)
%            + State-space trajectory selection: pod-greedy(0), pod(1), pca(2)
%            + Add x States per Iteration
%            + Orthogonalization Algorithm: orth(0), qr(1), svd(2)
%            + Solver: IRK3(0), Custom Solver with handle(-1)
%        (matrix,vector,scalar,handle)  ut
%        (vector)  x0 - initial state
%        (matrix)  yd - experimental data
%
% RETURNS:
%              (cell)  XP - {State-,Parameter-} Projection
%
% CITATION:
%    C. Himpe (2015). optmor - Optimization-Based Model Order Reduction (Version 2.1)
%    [Computer Software]. Available from http://???
%
% KEYWORDS:
%    model reduction, combined reduction, greedy
%*

    % Version Info
    if(nargin==1 && strcmp(f,'version')), XP = 2.1; return; end;

    % Default Arguments
    if(nargin<7  || isempty(co)), co = 1; end; % Assume unit covariance
    if(nargin<8  || isempty(nf)), nf = 0; end; % Assume default options
    if(nargin<9  || isempty(ut)), ut = 1; end; % Assume impulse input
    if(nargin<10 || isempty(x0)), x0 = 0; end; % Assume zero initial state
    if(nargin<11 || isempty(yd)), yd = 0; end; % Assume no experimental data

    % System Constants
    J = s(1); % Number of inputs
    N = s(2); % Number of states
    O = s(3); % Number of outputs

    h = t(2); % Time step width
    T = round((t(3)-t(1))/h); % Number of time steps
    r = abs(r); % Ensure positivity
    Q = numel(q); % Number of parameters

    % Discretize Procedural Input
    if(isa(ut,'function_handle'))
        uf = ut;
        ut = zeros(J,T);
        for l=1:T
            ut(:,l) = uf(l*h);
        end;
    end;

%% LAZY ARGUMENTS

    if(numel(co)==1), co = co*ones(Q,1); end
    if(numel(nf)==1), nf = zeros(1,12); end
    if(numel(ut)==1), ut = [ut*ones(J,1),sparse(J,T-1)]; end
    if(numel(x0)==1), x0 = x0*ones(N,1); end

%% SETUP

    % SET ABORT CRITERIA
    if(r >= 1)
        n = r;
    else
        n = N;
    end

    % SET DEFAULT TIKHONOV REGULARIZATION
    if(nf(3)==0),
        nf(3) = 0.5;
    end

    % CHECK FOR DATA
    if(nf(4) && ( (numel(yd)==1 && yd==0) || size(yd,1)~=O || size(yd,2)~=T) ),
        error('ERROR! optmor: yd data dimension mismatch!');
    end;

    % SET MINIMUM MONTE-CARLO BASE SIZE
    if(nf(5)),
        nf(5) = max(2,round(nf(5)));
    end

    % SEED RANDOMIZERS
    if(nf(6)),
        rand('seed',nf(6));
        randn('seed',nf(6));
    end;

    % SET MINIMUM STATE PRINCIPAL COMPONENTS
    if(nf(8)==0),
         nf(8) = 1;
    end

    % Precision Matrix
    if(size(co,2)==1),
        W = spdiags(1.0./co,0,Q,Q);
    else,
        W = pinv(full(co));
    end

    switch(nf(10)) % ODE Integrator

        case 0, % IRK3
            ode = @irk3;

        case -1, % Custom
            global CUSTOM_ODE;
            ode = CUSTOM_ODE;
    end

%% OBJECTIVE FUNCTIONAL SETUP

    b1 = nf(2);
    b2 = nf(3);
    bd = nf(4);

    if(b1==0 && b2==0 && bd==0), reg = 0; end; %  NONE

    if(b1~=0 && b2==0 && bd==0), reg = @(p,y) b1*norm1w(p,W); end % L1

    if(b1==0 && b2~=0 && bd==0), reg = @(p,y) b2*norm2w(p,W); end; % L2

    if(b1==0 && b2==0 && bd~=0), reg = @(p,y) bd*norm2t(yd-y,h); end; % DD

    if(b1~=0 && b2~=0 && bd==0), reg = @(p,y) b1*norm1w(p,W) + b2*norm2w(p,W); end; % L1 + L2

    if(b1~=0 && b2==0 && bd~=0), reg = @(p,y) b1*norm1w(p,W) + bd*norm2t(yd-y,h); end % L1 + DD

    if(b1==0 && b2~=0 && bd~=0), reg = @(p,y) b2*norm2w(p,W) + bd*norm2t(yd-y,h); end; % L2 + DD

    if(b1~=0 && b2~=0 && bd~=0), reg = @(p,y) b1*norm1w(p,W) + b2*norm2w(p,W) + bd*norm2t(yd-y,h); end; % L1 + L2 + DD

%% MAIN LOOP

    % Set Initial Parameter
    p = q;
    P = p;

    % Compute Trajectory for Initial Parameter
    z = ode(f,1,h,T,x0,ut,p);
    X = prn(z,nf(7),nf(8),1);

    opt = optimset('Display','off');

    for I=2:n

        fr = @(x,u,p) X'*f(X*x,u,p);
        gr = @(x,u,p) g(X*x,u,p);
        x0r = X'*x0;
        mmz = @(p,y) reg(p,y) - norm2t(y-ode(f,g,h,T,x0,ut,p),h);

        if(nf(5))
            mc = orth(rand(Q,nf(5)));
            p = ones(nf(5),1); % mc\q;
            J = @(p) mmz(mc*p,ode(fr,gr,h,T,x0r,ut,P*P'*mc*p));
        else
            p = q;
            J = @(p) mmz(p,ode(fr,gr,h,T,x0r,ut,P*P'*p));
        end

        switch(nf(1)) % Greedy Sampling Algorithm

            case 0, % Unconstrained (Quasi-Newton)
                p = fminunc(J,p,opt);

            case 1, % Derivative-Free (Nelder-Mead)
                p = fminsearch(J,p,opt);

            case -1, % Custom Optimizer
                global CUSTOM_OPT;
                p = CUSTOM_OPT(J,p);
        end

        if(nf(5)), p = mc*p; end;
        P = inc([P,p],nf(9));

        z = ode(f,1,h,T,x0,ut,p);
        z = z - X*(X'*z);
        x = prn(z,nf(7),nf(8),1);
        X = ([X,x]);

        if(mod(I,100)==0) fprintf('#'); elseif(mod(I,10)==0), fprintf('+'); else, fprintf('|'); end;

        if(r<1.0)
            yf = ode(f,g,h,T,x0,ut,q);
            yr = ode(@(x,u,p) X'*f(X*x,u,p),@(x,u,p) g(X*x,u,p),h,T,X'*x0,ut,P*(P'*q));
            if(norm2t(yf-yr,h)/norm2t(yf,h) < r), break; end;
        end
    end

    XP = {X,P};
end


%% ======== SQUARED NORMS ========
function x = norm2t(X,h)

    x = h*sum(X(:).*X(:));
end

function x = norm1w(X,w)

    x = sum(abs(w*X));
end

function x = norm2w(X,w)

    x = X'*X; %(X'*w')*(w*X);
end

function x = norm2(X)

     x = X'*X;
end

%% ======== STATE-SELECT ========
function v = prn(u,o,s,w)

    switch(o)
        case 1, % PCA
            [v,dummy,dummy2] = svds(bsxfun(@minus,u,mean(u,2),s));

        case 2, % Weighted POD
            [v,dummy,dummy2] = svds(w*u,s);

        otherwise, % = case 0, % POD
            [v,dummy,dummy2] = svds(u,s);
    end;
end


%% ======== INCORPORATE ========
function v = inc(u,o)

    switch(o)
        case 1, % QR
            [v,dummy] = qr(u,0);

        case 2, % SVD
            [v,dummy,dummy2] = svd(u,'econ');

        otherwise, % ORTH
            v = orth(u);
    end;
end

%% ======== DEFAULT ODE INTEGRATOR ========
function x = irk3(f,g,h,T,z,u,p)

    if(isnumeric(g) && g==1), g = @(x,u,p) x; end;

    k1 = h*f(z,u(:,1),p); % 2nd Order Midpoint RK2 for starting value
    k2 = h*f(z + 0.5*k1,u(:,1),p);
    z = z + k2;
    x(:,1) = g(z,u(:,1),p);

    x(end,T) = 0; % preallocate trajectory

    k1 = h*f(z,u(:,2),p); % 2nd start value (improves impulse response)
    k2 = h*f(z + 0.5*k1,u(:,2),p);
    z = z + k2;
    x(:,2) = g(z,u(:,2),p);

    for t=3:T % 3rd Order Improved Runge-Kutta IRK3
        l1 = h*f(z,u(:,t),p);
        l2 = h*f(z + 0.5*l1,u(:,t),p);
        z = z + (2.0/3.0)*l1 + (1.0/3.0)*k1 + (5.0/6.0)*(l2 - k2);
        x(:,t) = g(z,u(:,t),p);
        k1 = l1;
        k2 = l2;
    end;
end

