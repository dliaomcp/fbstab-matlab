% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% computes solutions of FBstab linear systems
% using a riccati recursion
classdef RiccatiLinearSolver < handle
	properties(Access = public)
	ztol = 100*eps;
	alpha = 0.95;
	data;
	sigma = sqrt(eps);

	mus;
	gamma;

	L;
	P;
	M
	SM;
	AM;
	Sigma;

	h;
	theta;
	end % properties

	methods(Access = public)

	function o = RiccatiLinearSolver(data)
		o.data = data;
		[nx,nu,nc,N] = data.OcpSize();

		o.P = zeros(nx,nu,N+1);
		o.Sigma = zeros(nu,nu,N+1);
		o.M = zeros(nx,nx,N+1);
		o.L = zeros(nx,nx,N+1);
		o.SM = zeros(nu,nx,N+1);
		o.AM = zeros(nx,nx,N+1);

		o.h = zeros(nx,N+1);
		o.theta = zeros(nx,N+1);
	end

	function Factor(o,x,xbar,sigma)
		[nx,nu,nc,N] = OcpSize(o.data);
		[nz,nl,nv] = ProblemSize(o.data);
		o.sigma = sigma;

		% Compute the barrier terms
		ys = x.y + sigma*(x.v - xbar.v);
		[gamma,mu] = o.dphi(ys,x.v);
		mus = sigma*gamma + mu;
		o.mus = mus;
		o.gamma = gamma;

		% Form the augmented Hessian
		Q = o.data.Q;
		R = o.data.R;
		S = o.data.S;
		B = gamma./mus;
		for i = 1:N+1
			Q(:,:,i) = Q(:,:,i) + sigma*eye(nx);
			R(:,:,i) = R(:,:,i) + sigma*eye(nu);

			Q(:,:,i) = Q(:,:,i) + ...
			o.data.E(:,:,i)' * o.diagmult(o.data.E(:,:,i),B(1+(i-1)*nc:i*nc));

			R(:,:,i) = R(:,:,i) + ...
			o.data.L(:,:,i)' * o.diagmult(o.data.L(:,:,i),B(1+(i-1)*nc:i*nc));

			S(:,:,i) = S(:,:,i) + ...
			o.data.L(:,:,i)' * o.diagmult(o.data.E(:,:,i),B(1+(i-1)*nc:i*nc));
		end

		% Temp variables
		Linv = zeros(nx,nx);

		% Begin the Riccati recursion
		% base case Pi = sigma I, L = chol(Pi)
		o.L(:,:,1) = sqrt(sigma)*eye(nx);
		for i = 1:N
			% Get the inverse factorization
			Linv = o.L(:,:,i)\eye(nx);

			% Compute QQ = Q + inv(Pi) and factor
			o.M(:,:,i) = Q(:,:,i) + Linv'*Linv;
			o.M(:,:,i) = chol(o.M(:,:,i),'lower');

			% Compute AM = A*inv(M)' and SM = S*inv(M)'
			o.AM(:,:,i) = o.data.Ak(:,:,i)/(o.M(:,:,i)');
			o.SM(:,:,i) = S(:,:,i)/(o.M(:,:,i)');

			% Compute Sigma = chol(R - S*inv(QQ)*S')
			o.Sigma(:,:,i) = ...
			R(:,:,i) - o.SM(:,:,i)*o.SM(:,:,i)';
			o.Sigma(:,:,i) = chol(o.Sigma(:,:,i),'lower');

			% Compute P = (A*inv(QQ)*S' - B)* inv(Sigma)';
			o.P(:,:,i) = o.AM(:,:,i)*o.SM(:,:,i)' - o.data.Bk(:,:,i);
			o.P(:,:,i) = o.P(:,:,i)/(o.Sigma(:,:,i)');

			% Compute Pi k+1
			o.L(:,:,i+1) = o.P(:,:,i)*o.P(:,:,i)' + ...
			o.AM(:,:,i)*o.AM(:,:,i)' + sigma*eye(nx);

			% Compute a cholesky factorization of Pi k+1
			o.L(:,:,i+1) = chol(o.L(:,:,i+1),'lower');

		end

		% Finish the recursion
		Linv = o.L(:,:,N+1)\eye(nx);
		% Compute QQ = Q + inv(Pi) and factor
		o.M(:,:,N+1) = Q(:,:,N+1) + Linv'*Linv;
		o.M(:,:,N+1) = chol(o.M(:,:,N+1),'lower');

		o.SM(:,:,N+1) = S(:,:,N+1)/(o.M(:,:,N+1)');

		o.Sigma(:,:,N+1) = ...
		R(:,:,N+1) - o.SM(:,:,N+1)*o.SM(:,:,N+1)';
		o.Sigma(:,:,N+1) = chol(o.Sigma(:,:,N+1),'lower');

	end % factor

	function Solve(o,r,dx)
		[nx,nu,nc,N] = OcpSize(o.data);
		[nz,nl,nv] = ProblemSize(o.data);

		% compute the reduced residuals
		r1 = r.rz - o.data.AT(r.rv./o.mus);
		r2 = -r.rl;

		LT.LT = true;
		UT.UT = true;
		% compute the forward recursion for h and theta
		o.theta(:,1) = r2(1:nx);
		% compute h
		o.h(:,1) = o.L(:,:,1)\o.theta(:,1);
		o.h(:,1) = o.L(:,:,1)'\o.h(:,1) - r1(1:nx);
		for i = 1:N
			ii = (i-1)*(nx+nu);
			iii = i*(nx+nu);
			ru = r1(ii+nx+1:ii+nx+nu);
			rxp = r1(1+iii:iii+nx);
			rlp = r2(1 + i*nx:(i+1)*nx);

			% compute theta(i+1)
			rt = ru+ o.SM(:,:,i)*linsolve(o.M(:,:,i),o.h(:,i),LT);
			rt = linsolve(o.Sigma(:,:,i),rt,LT);

			rt2 = linsolve(o.M(:,:,i),o.h(:,i),LT);

			o.theta(:,i+1) = o.P(:,:,i)*rt;
			o.theta(:,i+1) = o.theta(:,i+1) + o.AM(:,:,i)*rt2;
			o.theta(:,i+1) = o.theta(:,i+1) + rlp;

			% compute o.h(i+1)
			o.h(:,i+1) = linsolve(o.L(:,:,i+1),o.theta(:,i+1),LT);
			o.h(:,i+1) = linsolve(o.L(:,:,i+1)',o.h(:,i+1),UT) - rxp;
		end

		% terminal trio
		ii = (N)*(nx+nu);
		ru = r1(ii+nx+1:ii+nx+nu);

		uN = ru + o.SM(:,:,N+1)*(o.M(:,:,N+1)\o.h(:,N+1));
		uN = o.Sigma(:,:,N+1)\uN;
		uN = o.Sigma(:,:,N+1)'\uN;

		xN = o.M(:,:,N+1)\o.h(:,N+1) + o.SM(:,:,N+1)'*uN;
		xN = -o.M(:,:,N+1)'\xN;

		lN = o.theta(:,N+1) + xN;
		lN = o.L(:,:,N+1)\lN;
		lN = -o.L(:,:,N+1)'\lN;

		dx.z(ii+1:ii+nx+nu) = [xN;uN];
		dx.l(1+N*nx:(N+1)*nx) = lN;

		% back substituite for the full solution
		for i = N:-1:1
			ii = (i-1)*(nx+nu);
			ru = r1(ii+nx+1:ii+nx+nu);
			lp = dx.l(1+nx*i:nx*(i+1));

			uidx = ii+nx+1:ii+nx+nu;
			xidx = ii +1:ii+nx;
			lidx = 1+nx*(i-1):nx*(i);

			% compute u
			dx.z(uidx) = o.SM(:,:,i)* (o.M(:,:,i)\o.h(:,i));
			dx.z(uidx) = linsolve(o.Sigma(:,:,i),ru + dx.z(uidx),LT);
			dx.z(uidx) = dx.z(uidx) + o.P(:,:,i)'*lp;
			dx.z(uidx) = linsolve(o.Sigma(:,:,i)',dx.z(uidx),UT);

			% compute x
			dx.z(xidx) =  linsolve(o.M(:,:,i),o.h(:,i),LT);
			dx.z(xidx) = o.AM(:,:,i)'*lp + dx.z(xidx);
			dx.z(xidx) = o.SM(:,:,i)' * dx.z(uidx) + dx.z(xidx);
			dx.z(xidx) = -linsolve(o.M(:,:,i)',dx.z(xidx),UT);

			% compute lambda
			dx.l(lidx) = o.theta(:,i) + dx.z(xidx);
			dx.l(lidx) = linsolve(o.L(:,:,i),dx.l(lidx),LT);
			dx.l(lidx) = -linsolve(o.L(:,:,i)',dx.l(lidx),UT);
		end

		% recover ieq duals
		dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus;
		% dy for linesearch
		dx.y = o.data.b - o.data.A(dx.z);
	end


	function A = diagmult(o,A,d)
		[m,n] = size(A);
		for i = 1:m
			A(i,:) = d(i)*A(i,:);
		end
	end
	
	% compute an element of the C-differential 
	function [gamma,mu] = dphi(o,a,b,nv)
		[~,~,q] = ProblemSize(o.data);
		ztol = o.ztol;
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
				gamma(i) = o.alpha*(1- d);
				mu(i) = o.alpha*(1- d);
			elseif S1(i)
				gamma(i) = o.alpha*(1- a(i)/r(i)) + (1-o.alpha)*b(i);
				mu(i) = o.alpha*(1- b(i)/r(i)) + (1-o.alpha)*a(i);
			else % S2
				gamma(i) = o.alpha*(1- a(i)/r(i));
				mu(i) = o.alpha*(1- b(i)/r(i));
			end
		end
	end% dphi

	end % methods


end % riccati