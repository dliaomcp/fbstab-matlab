% fbstab_dense: Proximally Stabilized Semismooth Method for solving convex QPs
% [x,v,out] = fbstab_dense(qp,x0,v0,opts) attempts to solve the 
% following convex quadratic programming problem
%
% min.  1/2 x'Hx + f'x
% s.t.  Ax <= b
%
% Aside from convexity there are no assumptions made about the problem
% This method can detect unboundedness/infeasibility
% and can exploit arbitrary initial guesses. This implementation uses 
% dense linear algebra and dense factorizations.
%
% Inputs: 
% qp: A structure with the following fields
% 	H: n x n Symmetric positive semidefinite Hessian matrix
% 	f: n x 1 Forcing vector
% 	A: q x n Constraint Jacobian
% 	b: q x 1 Constraint vector
% x0: n x 1 Primal initial guess, use 0 if unsure
% v0: q x 1 Dual initial guess, use 0 if unsure
%
% opts: Options structure, see fbstab_options.m
%
% Outputs:
% x: n x 1 Primal solution
% v: q x 1 Dual solution
% If the QP is dual infeasible x will contain an unbounded 
% 	descent direction
% If the QP is primal infeasible v will contain a certificate of
% 	infeasibility
% out: Structure with the following fields
% 	prox_iters: number of proximal iterations taken
% 	newton_iters: total number of newton iterations
% 	res:  KKT residual
% 	eflag: Exit flag
% 		0: success
% 		-1: Maximum number of iterations exceeded
% 		-2: Problem is infeasible
% 		-3: Problem is unbounded below (dual infeasible)
% 		-4: Problem is primal and dual infeasible
%
% Details about the the FBstab algorithm can be found at:
% https://arxiv.org/abs/1901.04046
%
% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% Copyright 2018-2019 University of Michigan

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [x,v,out] = fbstab(qp,x0,v0,opts)

	if nargin < 4
		opts = struct();
	end

	n = length(qp.f);
	q = length(qp.b);

	% default parameters
	sigma = sqrt(eps); % initial stabilization parameter
	max_newton_iters = 100; % maximum number of allowable newton iterations
	max_prox_iters = 20;
	tol = 1e-6; % absolute error tolerance
	dtol = 1e-12; % stall tolerance
	rtol = 1e-12; % relative error tolerance
	enable_prescaling = false; % automatically rescale the QP
	inf_tol = 1e-8; % infeasibility tolerance
	check_infeasibility = true;
	veps_max = 1; % maximum error tolerance for the semismooth solver
	alpha = 0.95; % convex combination of FB and positive orthant penalization

	% overwrite defaults with options if they're given
	if isfield(opts,'sigma')
		sigma = opts.sigma;
	end
	if isfield(opts,'prox_iters')
		max_prox_iters = opts.prox_iters;
	end
	if isfield(opts,'tol')
		tol = opts.tol;
	end
	if isfield(opts,'rtol')
		rtol = opts.rtol;
	end
	if isfield(opts,'dtol')
		dtol = opts.dtol;
	end
	if isfield(opts,'enable_prescaling')
		enable_prescaling = opts.enable_prescaling;
	end
	if isfield(opts,'check_infeasibility')
		check_infeasibility = opts.check_infeasibility;
	end
	if isfield(opts,'alpha')
		alpha = opts.alpha;
	end
	if isfield(opts,'max_newton_iters')
		max_newton_iters = opts.max_newton_iters;
	end

	sc = struct();
	sc.pr = ones(n,1); % primal scaling
	sc.ieq = ones(q,1); % dual scaling
	sc.c = 1; % cost scaling

	E0 = nr(x0,v0,qp,sc,alpha);

	% precondition the problem
	if enable_prescaling
		[sc,qp] = ruiz_eq(qp,opts);
	end
	
	out = struct();
	out.eflag = -1;

	% initialize
	x = x0./sc.pr;
	v = v0./sc.ieq*sc.c;
	dx = ones(size(x0));
	dv = ones(size(v0));
	prox_iters = 0;
	newton_iters = 0;

	% pick the initial value of delta so that 
	% \eps_0 = min(E0,veps_max)
	delta = min(E0,veps_max)/sigma;

	iters = 0;
	for j = 1:max_prox_iters
		% convergence check
		E = nr(x,v,qp,sc,alpha);
		if E <= tol + E0*rtol || norm(dx)+norm(dv) <= dtol
			out.eflag = 0;
			break;
		elseif newton_iters > max_newton_iters
			break;
			out.eflag = -1;
		end
		
		% update inner solver tolerance
		delta = delta * (1/5);

		% safeguards
		inner_tol = min(E,delta*sigma);
		inner_tol = max(inner_tol,1e-13);

		iters = newton_iters;
		% approximately solve the proximal subproblem
		[xp,vp,R,newton_iters] = pfb(qp,x,v,x,v,sigma,inner_tol,alpha,newton_iters,max_newton_iters,opts);

		if iters - newton_iters== 0
			delta = delta/100;
		end

		%% Proximal refinement 
		% check if the proximal step makes sufficient progress, 
		% if not perform more newton iterations
		% it also exits if the tolerance is small enough
		% since if the problem is infeasible then reducing the
		% residual is not possible
		E = dist(x,v,qp,sc,alpha);
		for i = 1:5
			Ep = dist(xp,vp,qp,sc,alpha);
			if Ep < E || R <= min(100*eps,inner_tol)
				break;
			else
				% update inner solver tolerance
				delta = delta * (1/5);
				% safeguards
				inner_tol = min(E,delta*sigma);
				inner_tol = max(inner_tol,100*eps);
				% refine the solution using more Newton steps
				[xp,vp,R,newton_iters] = ...
				pfb(qp,xp,vp,x,v,sigma,inner_tol,alpha,newton_iters,max_newton_iters,opts);
			end
		end

		% compute increments and update solution
		dx = xp-x;
		dv = vp-v;
		x = xp;
		v = vp;

		% infeasibility detection
		if check_infeasibility
			% dual infeasibility
			w = infnorm(dx);
			c1 = infnorm(qp.H*dx) <= inf_tol * w;
			c2 = (qp.f)'*dx < 0;
			c3 = max(qp.A*dx) <= 0;
			dual_infeasible = c1 && c2 && c3 && (w >= eps);
			
			% primal infeasibility
			u = infnorm(dv);
			c4 = infnorm(qp.A'*dv) <= inf_tol*u;
			c5 = (qp.b)'*max(0,dv) < 0;
			primal_infeasible = c4 && c5 && (u >= eps);

			if dual_infeasible
				out.eflag = -3;
				x = dx;
				break;
			elseif primal_infeasible
				out.eflag = -2;
				v = dv;
				break;
			elseif dual_infeasible && primal_infeasible
				out.eflag = -4;
				x = dx;
				v = dv;
				break;
			end
		end

		% record diagnostics
		prox_iters = prox_iters +1;
	end

	out.prox_iters = prox_iters;
	out.newton_iters = newton_iters;
	out.res = E;
	out.solver = 'fbstab';

	% invert the scaling
	x = sc.pr.*x;
	v = 1/sc.c*sc.ieq.*v;

end


function E = dist(x,v,qp,sc,alpha)
	r1 = qp.H*x+qp.f + (qp.A)'*v;
	y = qp.b-qp.A*x;
	r2 = alpha*min(y,v) + (1-alpha)*max(0,y).*max(0,v);
	E = norm([r1;r2]);
end

function E = nr(x,v,qp,sc,alpha)
	r1 = qp.H*x+qp.f + (qp.A)'*v;
	y = qp.b-qp.A*x;
	r2 = min(y,v);
	E = norm([r1;r2]);
end


% solves the following GE
% Hx + f + A'v + sigma*(x-xbar) = 0
% b - Ax + sigma*(v-vbar) N_+(v) \ni 0
% using a semismooth penalized Fischer-Burmeister (pfb) method
function [x,v,E,niters] = pfb(qp,x0,v0,xbar,vbar,sigma,tol,alpha,niters,max_iters,opts)

	n = length(qp.f);
	q = length(qp.b);

	% parameters
	beta = 0.7; % backtracking parameter
	eta = 1e-8; % sufficient decrease parameter
	lnm = 5; % recurrence length for the non-monotone linesearch
	lsmax = 20; % maximum number of allowable linesearch iterations
	max_inner_iters = 200; % maximum number of Newton iterations allowed
	if isfield(opts,'eta')
		eta = opts.eta;
	end
	if isfield(opts,'lnm')
		lnm = opts.lnm;
	end
	if isfield(opts,'beta')
		beta = opts.beta;
	end
	if isfield(opts,'lsmax')
		lsmax = opts.lsmax;
	end
	if isfield(opts,'inner_iters')
		max_inner_iters = opts.inner_iters;
	end

	% initialization
	x = x0;
	v = v0;

	y = qp.b - qp.A*x + sigma*(v - vbar);
	% used to record past merit function 
	% values for the nonmonotone linesearch
	mrec = zeros(lnm,1); 

	% diagnostics
	E = 0;
	for j = 1:max_inner_iters
		% compute RHS
		r1 = -(qp.H*x + qp.f + qp.A'*v + sigma*(x - xbar));
		r2 = -phi(y,v,alpha);
		R = -[r1;r2];
        
        ER = norm(x-xbar) + norm(v-vbar);
		% convergence check
		E = norm(R);
		if E <= tol*min(1,ER) || niters >= max_iters
			break;
		end

		% compute iteration matrix
		Hs = qp.H + eye(n)*sigma;
		[gamma,mu] = dphi(y,v,alpha);
		mus = mu + sigma*gamma;
		% compute some diagonal matrix - matrix products efficiently
		B = zeros(size(qp.A));
		for i = 1:q
			B(i,:) = gamma(i)/mus(i)*qp.A(i,:);
		end

		% solve using reduced linear algebra
		% Computes: K = H + A'*diag(gamma./mus)*A;
		K = Hs + (qp.A)'*B;
		r = r1 - (qp.A)'*(r2./mus);

		LT.LT = true;
		UT.UT = true;

		K = chol(K,'lower');
		% dx = K\r;
		dx = linsolve(K,r,LT);
		% dx = (K')\dx;
		dx = linsolve(K',dx,UT);
		dv = (r2 + gamma.*(qp.A*dx))./mus;
		dy = -qp.A*dx + sigma*dv;

		dz = [dx;dv];

		% linesearch (non-monotone)
		% update record of merit function values
		mrec = circshift(mrec,1);
		% compute new merit function value
		mrec(1) = 1/2*E^2; 
		m0 = max(mrec); % largest recorded merit function value
		t = 1;% initialize stepsize
		for i = 1:lsmax
			xp = x+t*dx;
			yp = y+t*dy;
			vp = v+t*dv;

			Rp = ...
			[qp.H*xp + qp.f + qp.A'*vp + sigma*(xp - xbar);
			phi(yp,vp,alpha)];
			mp = 1/2*norm(Rp)^2;
			% this is the same as mp <= m0 + t*eta*(V'*R)'dz
			if mp <= m0 - 2*t*eta*mrec(1)
				break;
			else
				t = beta*t;
			end
		end
		% update
		x = x + t*dx;
		v = v + t*dv;
		y = y + t*dy;
		niters = niters + 1;
	end % end main loop

	% perform projection on the final output
	v = max(0,v);
end


% compute the complimentarity function
function y = phi(a,b,alpha)

	yfb = a + b - sqrt(a.^2 + b.^2);
	ypen = max(0,a).*max(0,b);

	y = alpha* yfb + (1-alpha)*ypen;

end

% compute an element of the C-differential 
function [gamma,mu] = dphi(a,b,alpha)
	q = length(a);
	ztol = 100*eps;
	% computes an element from the C differential
	r = sqrt(a.^2 + b.^2);
	S0 =  r <= ztol; % zero elements
	S1 = a > 0 & b> 0; % positive orthant
	S2 = ~(S1 | S0); % the rest

	gamma = zeros(q,1);
	mu = zeros(q,1);

	for i = 1:q
		if S0(i)
			d = 1/sqrt(2);
			gamma(i) = alpha*(1- d);
			mu(i) = alpha*(1- d);
		elseif S1(i)
			gamma(i) = alpha*(1- a(i)/r(i)) + (1-alpha)*b(i);
			mu(i) = alpha*(1- b(i)/r(i)) + (1-alpha)*a(i);
		else % S2
			gamma(i) = alpha*(1- a(i)/r(i));
			mu(i) = alpha*(1- b(i)/r(i));

		end
	end
end

function y = infnorm(v)
	y = max(abs(v));
end



% [sc,qp] = ruiz_equlibration(qp,opts)
% Compute a diagonal scaling and cost scaling c which 
% reduces the condition number of M = [H,A';A,0]
% The original problem is 
% min. 1/2 x'Hx + f'x
% s.t. Ax <= b
% initial conditions should be rescaled as
% x0 = x0./p, v0 = v0./d * c
% the solutions can be recovered as
% xopt = xopt.*p, vopt = 1/c * vopt.*d


function [sc,qp,out] = ruiz_eq(qp,opts)
	if nargin < 2
		opts = struct();
	end

	max_ruiz_iters = 30;
	ruiz_tol = 1e-2;
	scale_cost = false;

	if isfield(opts,'max_ruiz_iters')
		max_ruiz_iters = opts.max_ruiz_iters;
	end
	if isfield(opts,'ruiz_tol')
		ruiz_tol = opts.ruiz_tol;
	end
	if isfield(opts,'scale_cost')
		scale_cost = opts.scale_cost;
	end

	n = length(qp.f);
	q = length(qp.b);
	N = n+q;

	c = 1;
	delta = zeros(N,1);
	s = ones(N,1);
	gamma = 1;

	M = [qp.H,qp.A';qp.A,zeros(q)];
	temp = zeros(n,1);
	j = 0;
	while infnorm(1-delta) > ruiz_tol && j <= max_ruiz_iters
		for i = 1:N
			rownorm = infnorm(M(i,:));
			% safeguard against large scalings
			if rownorm < 1e-4
				rownorm = 1;
			end
			delta(i) = 1/sqrt(rownorm);
		end
		M = lrscale(M,delta,delta);
		qp.f = delta(1:n).*qp.f;
		qp.b = delta(n+1:end).*qp.b;

		if scale_cost
			for i = 1:n
				temp(i) = infnorm(M(i,1:n));
			end
			cnorm = max(mean(temp),infnorm(qp.f));
			% safeguard
			if cnorm < 1e-4
				cnorm = 1;
			end
			gamma = 1/cnorm;
			M(1:n,1:n) = gamma*M(1:n,1:n);
			qp.f = gamma*qp.f;
			c = gamma*c;
		end
		s = delta.*s;
		j = j+1;
	end
	qp.H = M(1:n,1:n);
	qp.A = M(n+1:end,1:n);
	p = s(1:n);
	d = s(n+1:end);

	sc.pr = p; % primal scaling
	sc.ieq = d; % dual scaling
	sc.c = c; % cost scaling
	out.iters = j;
end

% computes A = diag(l)*A*diag(r)
function A = lrscale(A,l,r)
	[m,n] = size(A);

	% scale the rows by l
	for i = 1:m
		A(i,:) = l(i)*A(i,:);
	end

	% scale the columns by r
	for i = 1:n
		A(:,i) = r(i)*A(:,i);
	end
end






