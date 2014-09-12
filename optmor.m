function XP = optmor(A,B,C,D,F,T,R,x,u,q,k,nf,yd)
% optmor (Version 1.1)
% by Christian Himpe, 2013-2014 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% About:
%  optmor (OPTimization-based Model Order Reduction)
%  reduces simultaneously parameters and states of
%  of a dynamic system of the following type:
%
%  x' = Ax + Bu + F
%  y  = Cx + Du
%
%  A in R^NxN
%  B in R^NxM
%  C in R^OxN
%  D in R^OxM
%  F in R^N
%
% Parameters:
%  (handle) A : returns system matrix: Ap = A(p) 
%  (matrix) B : input matrix (states x inputs)
%  (matrix) C : output matrix (outputs x states)
%  (matrix) D : feed-forward matrix (outputs x inputs)
%  (vector) F : source term (states x 1)
%  (vector) T : timing [start,step,stop]
%  (scalar) R : reduced state dimension or error bound
%  (vector) x : initial value (states x 1)
%  (matrix) u : input data (inputs x timesteps)
%  (matrix) k : prior covariance matrix (parameters x parameters)
%  (vector) nf : configuration (12 x 1) {optional}
% 		* reduction target (Dimension=0,Error=1) [not implemented yet]
% 		* state residual (POD-GREEDY=0,POD=1,POD_ORTH=2)
%		* optimization norm (INFTY=0,1=1,...)
%		* param orthogonalization type (QR=0,SVD=1,ORTH=2)
%		* number of state base additions
%		* Parameter Regularization Coefficient beta
%		* Data Regularization Parameter gamma
%		* use Data-Driven Regularization
%		* use Monte-Carlo Basis Enrichment
%		* Integration Order (currently only 4-5)
%		* Size of Monte-Carlo base
%		* check Parameter Stability
%  (matrix) yd : experimental data (outputs x timesteps) {optional}
%
% Output:
%  (cell)  XP : State, Parameter Projections (reduced x full)
%
% Keywords: Model Reduction, Combined Reduction, Bayesian Inversion
%
% Cite:
%	C. Himpe and M. Ohlberger;
%	"Data-Driven Combined State and Parameter Reduction for Inverse Problems";
%	Preprint arXiv.OC, : 2014
%
%*

global gt; % only for benchmarking pruposes
global gd; % only for benchmarking pruposes

% seed randomizers
randn('seed',1009);
rand('seed',1009);

% Determine System Dimensions
h = T(2);
t = (T(3)-T(1))/h;
N = sqrt(numel(A(q)));
J = size(B,2);
O = size(C,1);
Q = numel(q);


% Set Up Parameter Mapping
if(~isa(A,'function_handle')),

	N = sqrt(Q);
	A = @(p) reshape(p,[N N]);
end


% Expand Lazy Arguments
if(nargin<8),  x = 0; end;
if(nargin<9),  u = 1; end;
if(nargin<10), k = 1; end;
if(nargin<11), nf = [0,0,0,0,0,0,0,0,0,0,0,0]; end;
if(nargin<12), y = 0; end;

if(D==0), D = sparse(O,J); end;
if(numel(F)==1), F = F*ones(N,1); end;
if(numel(x)==1), x = x*ones(N,1); end;
if(numel(u)==1), u = u*[ones(J,1), sparse(J,t-1)]; end;
if(numel(k)==1), k = k*speye(numel(Q)); end;
if(numel(nf)<12), nf(1,12) = 0; end;


% Set Default Options
if(nf(5)==0),  nf(5)  = 1; end;
if(nf(6)==0),  bb = 0.1; else, bb = nf(6); end;
if(nf(7)==0),  cc = 0.1; else, cc = nf(7); end;
if(nf(11)==0), nf(11) = 2; end;


% Compute Precision Matrix
G = k;
G(k~=0) = 1.0./k(k~=0);


% Assemble Objective Function
norm1 = @(m) h*sum(sqrt(sum(m.*m)));
norm2 = @(m) h*sum(sum(m.*m));
norm8 = @(m) max(sqrt(sum(m.*m)));
normp = @(m) h*sum(sum(m.^nf(2)));

if(nf(8)==0)
	switch(nf(2))

		case 0,    j = @(p,y) bb*(p'*G*p) - norm8(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y);

		case 1,    j = @(p,y) bb*(p'*G*p) - norm1(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y);

		case 2,    j = @(p,y) bb*(p'*G*p) - norm2(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y);

		otherwise, j = @(p,y) bb*(p'*G*p) - normp(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y);
	end
else
	switch(nf(2))

		case 0,    j = @(p,y) bb*(p'*G*p) - norm8(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y) + cc*norm8(yd-y);

		case 1,    j = @(p,y) bb*(p'*G*p) - norm1(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y) + cc*norm1(yd-y);

		case 2,    j = @(p,y) bb*(p'*G*p) - norm2(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y) + cc*norm2(yd-y);

		otherwise, j = @(p,y) bb*(p'*G*p) - normp(ode(A(p),B,C,D,F,h,t,x,u,nf(10))-y) + cc*normp(yd-y);
	end		
end;


% Init Projections
X = prd(ode(A(q),B,speye(N),sparse(N,J),F,h,t,x,u,nf(10)),nf(5));
P = 0.1*randn(Q,1)+q;
p = q;

XP = {X,P};

gd(1,1,1) = toc(gt);	 % only for benchmarking pruposes


% Main Loop
for I=2:R,
	br = X'*B;
	fr = X'*F;
	cr = C*X;
	xr = X'*x;

	if(nf(9)==1),
		% Monte-Carlo Basis Enrichment
		qq = 0.1*randn(Q,nf(11))+q*ones(1,nf(11));

		PPqq = P*(P'*qq);
		kk = fminunc(@(k) j(qq*k,ode(X'*A(PPqq*k)*X,br,cr,D,fr,h,t,xr,u,nf(10))),ones(nf(11),1)./nf(11) ,optimset('Display','off'));
		p  = qq*kk;
	else,
		% Original Method
		p = fminunc(@(k) j(k,ode(X'*A(P*(P'*k))*X,br,cr,D,fr,h,t,xr,u,nf(10))),q ,optimset('Display','off'));
	end;

	% Test Stability
	if(eigs(A(p),1,'lr')>=0 && nf(12)), fprintf('8'); end

	% Append State Projections
	e = ode(A(p),B,speye(N),sparse(N,J),F,h,t,x,u,nf(10));
	f = ode(X'*A(p)*X,br,X,sparse(N,J),fr,h,t,xr,u,nf(10));

	switch(nf(3)),
		case 0,
			X = [X,prd(e-f,nf(5))];
		case 1,
			X = [X,prd(e,nf(5))];
		case 2,
			P = orth([X,prd(e,nf(5))]);
	end;

	% Append Parameter Projection
	switch(nf(4)),
		case 0,
			[P,dummy] = qr([P,p],0);
		case 1,
			[P,dummy,dummy2] = svd([P,p],'econ');
		case 2,
			P = orth([P,p]);
	end;

	XP = [XP;{X,P}];

	gd(1,1,I) = toc(gt);
end;

end


%%%%%%%% PRINCIPAL DIRECTION %%%%%%%%
function u = prd(x,n)

	[u,D,V] = svds(x,n);
end

%%%%%%%% ODE SOLVER %%%%%%%%
function y = ode(A,B,C,D,F,h,L,x,u,O)

	H = 1.0./h;

	if(exist('OCTAVE_VERSION')),
		f = @(y,t) A*y + B*u(:,1.0+min(round(t*H),L-1)) + F;
		y = C*lsode(f,x,linspace(0,h*L,L))' + D*U;
	else,
		f = @(t,y) A*y + B*u(:,1.0+min(round(t*H),L-1)) + F;
		y = C*deval( ode45(f,[0,h*L],x), linspace(0,h*L,L) ) + D*u;
	end;
end

