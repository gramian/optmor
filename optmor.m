function XP = optmor(f,g,s,t,r,q,nf,ut,x0,co,yd)
% optmor (Version 2.5)
% by Christian Himpe, 2013-2016 ( wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNTAX:
%    W = optmor(f,g,s,t,r,q,[nf],[ut],[x0],[co],[yd]);
%
% SUMMARY:
%    optmor - optimization-based model order reduction,
%    for the computation of combined state and parameter
%    reduced order models of state-space input-output systems.
%    Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%   (func handle)  f - system function handle; signature: xdot = f(x,u,p)
%   (func handle)  g - output function handle; signature:    y = g(x,u,p)
%        (vector)  s - system dimensions [inputs,states,outputs]
%        (vector)  t - time discretization [step,stop]
%        (scalar)  r - reduced order(>1) or error threshold(<1)
%        (vector)  q - nominal parameter
%        (vector,scalar) [nf = 0] - options, 6 components:
%            + Optimization Algorithm: fminunc(0), fminsearch(1), custom(-1)
%            + Lasso Regularization Weight: default(0)
%            + Tikhonov Regularization Weight: default(0.1)
%            + Data-Driven Regularization Weight: default(0)
%            + Number of Maximum Optimizer Iterations: default(4)
%            + Initial Parameter: last(0), random(1)
%        (matrix,vector,scalar,handle)  ut - input; default: delta impulse
%        (vector,scalar)  x0 - initial state; default: zeros
%        (matrix,vector,scalar)  [co = 1] - covariance matrix: unit
%        (matrix)  yd - experimental data: empty
%
% RETURNS:
%              (cell)  XP - {State-,Parameter-} Projection
%
% CITATION:
%    C. Himpe (2016). optmor - Optimization-Based Model Order Reduction
%    (Version 2.5) [Software]. Available from http://gramian.github.io/optmor .
%    doi: 10.5281/zenodo.46683 .
%
% KEYWORDS:
%    model reduction, combined reduction, greedy sampling
%*

    % Custom Solver
    global ODE;
    if(isa(ODE,'function handle')==0), ODE = @rk2; end;

    % Version Info
    if( nargin==1 && strcmp(f,'version') ), XP = 2.5; return; end;

    % Test if Octave
    if(~exist('OCTAVE_VERSION','builtin')), vec = @(m) m(:); end;

    % Default Arguments
    if( nargin<7  || isempty(nf) ), nf = 0; end; % Assume default options
    if( nargin<8  || isempty(ut) ), ut = 1.0; end; % Assume impulse input
    if( nargin<9  || isempty(x0) ), x0 = 0; end; % Assume zero initial state
    if( nargin<10 || isempty(co) ), co = 1.0; end; % Assume covariance
    if( nargin<11 || isempty(yd) ), yd = 0; end; % Assume no experimental data

    % System Constants
    J = s(1);               % Number of inputs
    N = s(2);               % Number of states
    O = s(3);               % Number of outputs
    h = t(1);               % Time step width
    T = floor(t(2)/h) + 1;  % Number of time steps
    Q = numel(q);           % Number of parameters

    % Linear Chirp Input
    if( isnumeric(ut) && numel(ut)==1 && ut==Inf )
        ut = @(t) 0.5*cos(pi*(t+10*t.*t))+0.5;
    end;

    % Discretize Procedural Input
    if(isa(ut,'function_handle'))
        uf = ut;
        ut = zeros(J,T);
        for l=1:T
            ut(:,l) = uf(l*h);
        end;
    end;

    % Lazy arguments
    if(numel(nf)<6),  nf(6) = 0; end;
    if(numel(ut)==1), ut(1:J,1) = ut./h; end;
    if(numel(x0)==1), x0(1:N,1) = x0; end;
    if(numel(co)==1), co(1:Q,1) = co; end;

    if(size(ut,2)==1), ut(:,2:T) = 0.0; end;
    if(size(co,1)==1 || size(co,2)==1), co = spdiags(co,0,Q,Q); end;

%% SETUP

    % Set Abort Critera
    r = abs(r);
    if(r >= 1.0)
        n = r;
    else
        n = Q;
    end;

    % Set Default Tikhonov Regularization Coefficient
    if(nf(3)==0 && nf(2)==0)
        nf(3) = 0.1;
    end;

    % Check for Data if Data-Driven Regularization
    if(nf(4) && ( (numel(yd)==1 && yd==0) || size(yd,1)~=O || size(yd,2)~=T) )
        error('ERROR! optmor: yd data dimension mismatch!');
    end;

    % Set Default Maximum Optimizer Iterations
    if(nf(5)==0)
        nf(5) = 4;
    end;

    % Set up Regularization Operators
    if(nf(2)~=0 && nf(3)==0 && nf(4)==0)
        R = @(p,y) nf(2)*norm(p,1);
    elseif(nf(2)==0 && nf(3)~=0 && nf(4)==0)
        R = @(p,y) nf(3)*norm(p,2)^2;
    elseif(nf(2)~=0 && nf(3)~=0 && nf(4)==0)
        R = @(p,y) nf(2)*norm(p,1) + nf(3)*norm(p,2)^2;
    elseif(nf(2)==0 && nf(3)==0 && nf(4)~=0)
        R = @(p,y) nf(4)*t(1)*norm(vec(y-yd),2)^2;
    elseif(nf(2)~=0 && nf(3)==0 && nf(4)~=0)
        R = @(p,y) nf(2)*norm(p,1) + nf(4)*t(1)*norm(vec(y-yd),2)^2;
    elseif(nf(2)==0 && nf(3)~=0 && nf(4)~=0)
        R = @(p,y) nf(3)*norm(p,2)^2 + nf(4)*t(1)*norm(vec(y-yd),2)^2;
    elseif(nf(2)~=0 && nf(3)~=0 && nf(4)~=0)
        R = @(p,y) nf(2)*norm(p,1) + nf(3)*norm(p,2)^2 ...
                                   + nf(4)*t(1)*norm(vec(y-yd),2)^2;
    end;

%% INIT LOOP

    fprintf('optmor progress:\n');

    % Set Initial Parameter Projection
    p = q;
    P = q./norm(q,2);

    % Compute Trajectory for Initial Parameter
    z = ODE(f,1,t,x0,ut,p);
    [X,dtemp,vtemp] = svds(z,1);

    % Set Optimizer Options
    flags = optimset('Display','off','MaxIter',nf(5));

    % Per Interation Timing
    global TIMINGS;
    TIMINGS = zeros(Q,1);

    fprintf('|');

%% MAIN LOOP

    for I=2:n

        k = tic;

        % Current reduced order model
        fr = @(x,u,p) X'*f(X*x,u,p);
        gr = @(x,u,p) g(X*x,u,p);
        x0r = X'*x0;

        % Termination Test
        if(r < 1.0 && t(1)*norm(vec(ODE(f,g,t,x0,ut,p) ...
                                   -ODE(fr,gr,t,x0r,ut,P*(P'*p))),2) < r)
            break;
        end;

        % Set up cost function
        jf = @(p,y) -t(1)*norm(vec(y-ODE(fr,gr,t,x0r,ut,P*(P'*p))),2)^2 +R(p,y);
        Jf = @(p) jf(p,ODE(f,g,t,x0,ut,p));       

        % Initial Parameter
        if(nf(6))
            p = q + co * randn(Q,1);
        end;

        % Greedy Sampling Algorithm
        switch(nf(1))

            case 0, % Unconstrained (Quasi-Newton)
                p = fminunc(Jf,p,flags);

            case 1, % Derivative-Free (Nelder-Mead)
                p = fminsearch(Jf,p,flags);

            case -1, % Custom Optimizer
                global FMIN;
                p = FMIN(Jf,p);
        end;

        % Extract and Incorporate New State Base
        z = ODE(f,1,t,x0,ut,p);
        z = z - X * (z' * X)';
        [x,dtemp,vtemp] = svds(z,1);
        X = gramschmidt(X,x);

        % Incorporate New Parameter Base
        P = gramschmidt(P,p);

        TIMINGS(I) = toc(k);

        if(mod(I,100)), fprintf('|'); else, fprintf('#'); end;
    end;

    fprintf('\n');

    XP = {X,P};
end

%% ======== RE-ITERATED GRAM-SCHMIDT ========
function Q = gramschmidt(Q,v)

    m = size(v,2);

    for I=1:m
        if(size(Q,2)>=size(Q,1))
            return;
        end;

        w = v(:,I);
        b = 0;

        for J=1:10

            w = w - Q * (w' * Q)';
            b = norm(w,2);
            w = w./b;
            if(b>0.1)
                break;
            end;
        end;

        Q = [Q,w];
    end;
end

%% ======== DEFAULT ODE INTEGRATOR ========
function x = rk2(f,g,t,z,u,p)

    if(isnumeric(g) && g==1), g = @(x,u,p) x; end;

    h = t(1);
    L = floor(t(2)/h) + 1;

    x(:,1) = g(z,u(:,end),p);
    x(end,L) = 0; % preallocate trajectory

    for l=2:L % 2nd order Ralston's Runge-Kutta Method
        k1 = h*f(z,u(:,l-1),p);
        k2 = h*f(z + 0.666666666666667*k1,u(:,l-1),p);
        z = z + 0.25*k1 + 0.75*k2;
        x(:,l) = g(z,u(:,l-1),p);
    end;
end

