% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% computes solutions of FBstab linear systems
% using a banded ldlt factorization
classdef ldlt_ms < handle
	
	properties(Access = public)

	ztol = 100*eps;
	alpha = 0.95;
	data; % pointer to the problem data

	% private later
	K; % KKT matrix
	mus; % barrier terms
	gamma;% barrier term
	sigma; 
	p;
	pinv;

	end % properties
	
	methods(Access = public)
	% constructor
	function o = ldlt_ms(data)
		o.data = data;
		[nz,nl,nv] = data.opt_sz();
		o.sigma = 0;
		o.mus = zeros(nv,1);
		o.gamma = zeros(nv,1);
		o.p = zeros(nl+nv,1);
		o.pinv = zeros(nl+nv,1);
		o.K = zeros(nl+nv,nl+nv);
	end

	% form and factor the KKT matrix
	function factor(o,x,xbar,sigma)
		[nx,nu,nc,N] = sz(o.data);
		[nz,nl,nv] = opt_sz(o.data);

		o.sigma = sigma;
		% compute the barrier terms
		ys = x.y + sigma*(x.v - xbar.v);
		[gamma,mu] = o.dphi(ys,x.v);
		mus = sigma*gamma + mu;
		o.mus = mus;
		o.gamma = gamma;

		% form the augmented Hessian
		Q = o.data.Q;
		R = o.data.R;
		S = o.data.S;
		B = gamma./mus;
		for i = 1:N+1
			% compute the barrier terms
			Q(:,:,i) = Q(:,:,i) + sigma*eye(nx);
			R(:,:,i) = R(:,:,i) + sigma*eye(nu);
			for j = 1:nc
				k = j + (i-1)*nc;
				% outer products
				Q(:,:,i) = Q(:,:,i) + ...
				B(k)*o.data.E(j,:,i)'*o.data.E(j,:,i);

				S(:,:,i) = S(:,:,i) + ...
				B(k)*o.data.L(j,:,i)'*o.data.E(j,:,i);

				R(:,:,i) = R(:,:,i) + ...
				B(k)*o.data.L(j,:,i)'*o.data.L(j,:,i);
			end
		end
		
		% form the KKT matrix
		o.K = diag([zeros(nz,1);-sigma*ones(nl,1)]);

		for i = 1:N+1
			im = (i-1)*(nx+nu) + 1;
			ip = i*(nx+nu);

			o.K(im:ip,im:ip) = ...
			[Q(:,:,i),S(:,:,i)';
			S(:,:,i),R(:,:,i)];
		end

		% G matrix
		o.K(nz+1:nz+nx,1:nx) = -eye(nx);
		o.K(1:nx,nz+1:nz+nx) = -eye(nx);
		for i = 2:N+1
			im = (i-1)*nx+1 + nz;
			ip = i*nx + nz;

			jm = (i-2)*(nx+nu)+1;
			jp = jm + 2*nx+nu-1;

			o.K(im:ip,jm:jp) = ...
			[o.data.Ak(:,:,i-1),o.data.Bk(:,:,i-1),-eye(nx)];

			o.K(jm:jp,im:ip) = o.K(im:ip,jm:jp)';
		end
		% permute
		[o.p,o.pinv] = o.ILperm();
		% factor the KKT matrix
		o.K = o.ldlt(o.K,o.p);
	end % factor

	function solve(o,r,dx)
		[nx,nu,nc,N] = sz(o.data);
		[nz,nl,nv] = opt_sz(o.data);

		% compute the reduced residuals
		r1 = r.rz - o.data.AT(r.rv./o.mus);
		r2 = -r.rl;
		R = [r1;r2];

		% permute
		p = o.p;
		% Solving using the factored KKT matrix
		% PKP' = LDL'
		n = nz+nl;
		bw = 2*nx+nu-1;
		% forward solve
		% R <- L\R;
		for j = 1:n
			for i = j+1:min(j+bw,n)
				R(p(i)) = R(p(i)) - o.K(p(i),p(j))*R(p(j));
			end
		end
		% diagonal solve
		% R <- D\R;
		for i = 1:n
			R(p(i)) = R(p(i))/o.K(p(i),p(i));
		end
		% backsolve
		% R <- L'\R;
		for j = n:-1:1
			for i = max(1,j-bw):j-1
				R(p(i)) =  R(p(i)) - o.K(p(j),p(i))*R(p(j));
			end
		end

		dx.z = R(1:nz);
		dx.l = R(nz+1:n);

		% recover ieq duals
		dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus;
		% dy for linesearch
		dx.y = o.data.b -o.data.A(dx.z) + o.sigma*dx.v;
	end % solve

	% banded ldlt factorization
	function A = ldlt(o,A,p)
		[nx,nu,nc,N] = o.data.sz();

		n = size(A,1);
		v = zeros(n,1);
		pp = zeros(n,1);
		t = 0;

		bw = 2*nx+nu-1;
		for k = 1:n
			% diagonal element
			t = A(p(k),p(k));
			ti = 1/t;
			if abs(t) < o.ztol
				disp('Matrix is near singular....');
			end
			m = min(k+bw,n);
			for i = k+1:m
				% save the column below the diagonal
				v(i) = A(p(i),p(k));
				A(p(i),p(k)) = v(i)*ti;
			end

			% propagate the factorizaton
			for i = k+1:m
				for j = k+1:i
					A(p(i),p(j)) = A(p(i),p(j)) - v(i)*v(j)*ti;
				end
			end
		end % main loop
	end % ldlt

	% computes the interleaving permutation
	% so that e.g., [x1;u1;x2;u2;l1;l2] -> [l1;x1;u1;l2;x2;u2]
	function [pp,pinv] = ILperm(o)
		[nx,nu,nc,N] = o.data.sz();
		nz = (nx+nu)*(N+1);
		nl = nx*(N+1);

		p = 1:nz+nl;
		pp = zeros(nz+nl,1);
		pinv = zeros(nz+nl,1);

		for i = 1:N+1
			for j = 1:nx
				% lambda
				k = (i-1)*(2*nx+nu) + j;
				pp(k) = p(nz+j + (i-1)*nx);

				% x
				k = nx+ (i-1)*(2*nx+nu) + j;
				pp(k) = p(j +(i-1)*(nx+nu));
			end
				% u
			for j = 1:nu
				k = 2*nx+ (i-1)*(2*nx+nu) + j;
				pp(k) = p(nx+j + (i-1)*(nx+nu));
			end
		end
		% inverse permutation
		for i = 1:length(pp)
			pinv(pp(i)) = i;
		end
	end
	% compute an element of the C-differential 
	function [gamma,mu] = dphi(o,a,b)
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
	end

	end % methods



end