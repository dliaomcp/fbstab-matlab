% fbstab_mpc: Proximally Stabilized Semismooth Method for solving QPs
% [x,out] = fbstab_mpc(x0,mpc,opts) attempts to solve the 
% following quadratic programming problem:
%
% min.  \sum_{i=0}^N 1/2 [x(i)]' * [Q(i) S(i)'] [x(i)] + [q(i)]'*[x(i)]
%                        [u(i)]    [S(i) R(i) ] [u(i)]   [r(i)]  [u(i)]
% s.t.  x(i+1) = A(i)*x(i) + B(i) u(i) + c(i), i = 0 ... N-1
%       x(0) = x0
%       E(i)*x(i) + L(i)*u(i) + d(i) <= 0,     i = 0 ... N
%
% Where [Q(i),S(i)';S(i),R(i)] is positive semidefinite
% 
% Aside from convexity there are no assumptions made about the problem
% This method can detect unboundedness/infeasibility
% and exploit arbitrary initial guesses. 
%
% Inputs: 
% x0:  A structure with the following fields
%	z: primal initial guess, arranged as z = [x0;u0,...,uN;xN], 
%	l: co-state initial guess
%	v: inequality dual initial guess
%
% mpc: A structure with the following fields
% 	Q = Q0,Q1,... \in [nx,nx,N+1]
% 	R = R0,R1,... \in [nu,nu,N+1]
% 	S = S0,S1,... \in [nu,nx,N+1]
%   q = q0,q1,... \in [nx,N+1]
%	r = r0,r1,... \in [nu,N+1]
%	
%	A = A0,A1,... \in [nx,nx,N]
%	B = B0,B1,... \in [nx,nu,N]
%   c = c0,c1,... \in [nx,N]
%	x0 \in [nx,1]
%
%	E = E0,E1,... \in [nc,nx,N+1]
%	L = L0,L1,... \in [nc,nu,N+1]
%	d = d0,d1,... \in [nc,N+1]
% 
% opts: A structure containing any of the following fields
% 	sigma{sqrt(eps)}: Initial stabilization parameter
%	max_newton_iters{500}: Maximum total number of Newton iterations
% 	max_prox_iters{100}: maximum number of prox iterations allowed
% 	max_inner_iters{100}: Maximum allowable number of inner iterations
% 	tol{1e-6}: Absolute tolerance
%   rtol{1e-12}: Relative tolerance
%   inftol{1e-8}: Infeasibility tolerance
%	dtol{1e-12}: Stall tolerance
% 	beta{0.7}: backtracking linesearch parameter 
% 	eta{1e-8}: sufficient decrease parameter
% 	lsmax{20}: maximum number of linesearch iterations
% 	alpha{0.95}: penalized FB function parameter
%	check_infeasibility{true}: check for feasibility
%	itol_max{1e-1}: maximum inner tolerance
%	itol_min{10*eps}: minimum inner tolerance
%	itol_red_factor{1/10}: reduction factor for the inner tolerance
%
%	display_level{0}: controls printed output
%				   0: Silent
%				   1: Final 
%				   2: Iter 
%				   3: Iter detailed 
%
%	linear_solver{ric}: Which linear solver to use
%	               ric: Ricatti recursion based method
%				   pcg: Implicit conjugate gradient
%				   ldl: Banded LDL
% 
% Outputs:
% x: Structure with the following fields
%	z: primal solution z = [x0;u0;...;uN;xN]
% 	l: costates
% 	v: inequality duals
% 
% u: The MPC control action (u0)
% 
% out: Structure with the following fields
% 	prox_iters: number of proximal iterations 
% 	newton_iters: number of newton iterations
% 	res:  KKT residual
% 	eflag: Exit flag
% 		 0: success
% 		-1: Maximum number of iterations exceeded
%		-2: Algorithm stalled
% 		-3: Problem is infeasible
% 		-4: Problem is unbounded below (dual infeasible)
% 		-5: Problem is primal and dual infeasible
%
% Details about the the FBstab algorithm can be found at:
% https://arxiv.org/abs/1901.04046
%
% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

%Copyright 2018-2019 University of Michigan

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [x,u,out] = fbstab_mpc(x0,mpc,opts)
	if nargin < 3
		opts = struct();
	end

	% default linear solver is the riccati method
	algo = 'ric';
	if isfield(opts,'linear_solver')
		algo = opts.linear_solver;
	end

	switch algo
		case 'ric'
			[x,out] = fbstab_ricatti(x0,mpc,opts);
		case 'pcg'
			[x,out] = fbstab_pcg(x0,mpc,opts);
		case 'ldl'
			[x,out] = fbstab_ldl(x0,mpc,opts);
		otherwise
			[x,out] = fbstab_ricatti(x0,mpc,opts);
	end
	[nx,nu] = size(mpc.B(:,:,1));
	u = x.z(nx+1:nx+nu);

end


function [x,out] = fbstab_ricatti(x0,mpc,opts)

	data = data_ms(mpc.Q,mpc.R,mpc.S,mpc.q,mpc.r,mpc.A,...
		mpc.B,mpc.c,mpc.x0,mpc.E,mpc.L,mpc.d);

	linsys = ricatti(data);
	feas = feas_ms(data);
	r1 = res_ms(data);
	r2 = res_ms(data);

	x1 = var_ms(data);
	x2 = var_ms(data);
	x3 = var_ms(data);
	x4 = var_ms(data);
	% solver object
	fbstab = fbstab_algo(data,linsys,feas,r1,r2,x1,x2,x3,x4,opts);

	[x,out] = fbstab.solve(x0);
end


function [x,out] = fbstab_pcg(x0,mpc,opts)
	% create QP data structure
	data = data_ss(mpc.Q,mpc.R,mpc.S,mpc.q,mpc.r,mpc.A,...
		mpc.B,mpc.c,mpc.x0,mpc.E,mpc.L,mpc.d);

	% get sizes
	[nx,nu,nc,N] = data.sz();

	% components objects
	linsys = pcg_ss(data);
	feas = feas_ss(data);

	r1 = res_ss(data);
	r2 = res_ss(data);

	x1 = var_ss(data);
	x2 = var_ss(data);
	x3 = var_ss(data);
	x4 = var_ss(data);

	% solver object
	fbstab = fbstab_algo(data,linsys,feas,r1,r2,x1,x2,x3,x4,opts);

	% extract the controls
	z = reshape(x0.z,[nx+nu,N+1]);
	u = z(nx+1:nx+nu,:);

	% input for the condensed solver
	y0 = struct();
	y0.v = x0.v;
	y0.l = x0.l;
	y0.z = u(:);

	% call the solver
	[y,out] = fbstab.solve(y0);
	u = y.z;
	uu = reshape(u,[nu,N+1]);

	x = struct();
	
	x.z = zeros((nx+nu)*(N+1),1);
	x.l = zeros(nx*(N+1),1);
	x.v = y.v;

	% compute the state sequence
	xx = recover_states(uu,data);
	nz = (nx+nu)*(N+1);
	x.z = reshape([xx;uu],[nz,1]);

	% recover the costates
	x.l = recover_costates(x.z,x.v,data);
end

% recover the states after solving the condensed problem
function x = recover_states(u,data)
	[nx,nu,nc,N] = data.sz();

	x = zeros(nx,N+1);
	x(:,1) = data.xt;
	for i = 1:N
		x(:,i+1) = data.Ak(:,:,i)*x(:,i) + data.Bk(:,:,i)*u(:,i) ...
		+ data.ck(:,i);
	end
end

% recover the co-states (eq duals)
% after solving the condensed problems
function l = recover_costates(z,v,data)
	[nx,nu,nc,N] = data.sz();

	% compute the lagrangian of the multiple shooting problem
	pp = -(data.f_ms() + data.H_ms(z) + data.AT_ms(v));
	p = reshape(pp,[nx+nu,N+1]);

	l = zeros(nx,N+1);
	l(:,N+1) = -p(1:nx,N+1);

	for i = N:-1:1
		l(:,i) = data.Ak(:,:,i)'*l(:,i+1) - p(1:nx,i);
	end
	l = l(:);
end

function [x,out] = fbstab_ldl(x0,mpc,opts)

	data = data_ms(mpc.Q,mpc.R,mpc.S,mpc.q,mpc.r,mpc.A,...
		mpc.B,mpc.c,mpc.x0,mpc.E,mpc.L,mpc.d);

	linsys = ldlt_ms(data);
	feas = feas_ms(data);
	r1 = res_ms(data);
	r2 = res_ms(data);

	x1 = var_ms(data);
	x2 = var_ms(data);
	x3 = var_ms(data);
	x4 = var_ms(data);
	% solver object
	fbstab = fbstab_algo(data,linsys,feas,r1,r2,x1,x2,x3,x4,opts);

	% call the solver
	[x,out] = fbstab.solve(x0);

end