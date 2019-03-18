% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% computes solutions of FBstab linear systems
% using backslash
classdef backslash < handle

properties(Access = public)
	
	ztol = 100*eps;
	K;

	alpha;
	sigma; 
	data; % pointer to the problem data
	
end% properties

methods(Access = public)
	% constructor

	function o = backslash(data)
		o.data = data;
		o.sigma =0;
	end

	% form and factor the KKT matrix
	% for this implementation you just form the KKT matrix
	function factor(o,x,xbar,sigma)
		[nx,nu,nc,N] = sz(o.data);
		[nz,nl,nv] = opt_sz(o.data);
		o.sigma = sigma;

		H = zeros(nz,nz);
		for i = 1:N+1
			im = (i-1)*(nx+nu) + 1;
			ip = i*(nx+nu);
			H(im:ip,im:ip) = ...
			[o.data.Q(:,:,i),o.data.S(:,:,i)';
			o.data.S(:,:,i),o.data.R(:,:,i)];
		end
		H = H + sigma*eye(nz);

		% G matrix
		G = zeros(nl,nz);
		G(1:nx,1:nx) = -eye(nx);
		for i = 2:N+1
			im = (i-1)*nx+1;
			ip = i*nx;

			jm = (i-2)*(nx+nu)+1;
			jp = jm + 2*nx+nu-1;

			G(im:ip,jm:jp) = ...
			[o.data.Ak(:,:,i-1),o.data.Bk(:,:,i-1),-eye(nx)];
		end

		% A matrix
		A = zeros(nv,nz);
		for i = 1:N+1
			im = 1+(i-1)*nc;
			ip = i*nc;

			jm = 1+(i-1)*(nx+nu);
			jp = i*(nx+nu);

			A(im:ip,jm:jp) = [o.data.E(:,:,i),o.data.L(:,:,i)];
		end

		% compute FB products
		ys = x.y + sigma*(x.v - xbar.v);
		[gamma,mu] = o.dphi(ys,x.v);
		mus = sigma*gamma + mu;

		CA = zeros(size(A));
		for i = 1:nv
			CA(i,:) = gamma(i)*A(i,:);
		end

		% form the KKT matrix
		o.K = ...
		[H,G',A';
		-G,sigma*eye(nl),zeros(nl,nv);
		-CA,zeros(nv,nl),diag(mus)];
	end

	% solve the equation
	% V*dx = r
	% r is a residual object
	% dx is the search direction
	function solve(o,r,dx)
		R = [r.rz;r.rl;r.rv];

		x = o.K\R;

		[nz,nl,nv] = opt_sz(o.data);
		dx.z = x(1:nz);
		dx.l = x(nz+1:nz+nl);
		dx.v = x(nz+nl+1:end);

		% compute an appropriate dy for the linesearch
		dx.y = o.data.b -o.data.A(dx.z);

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


end % 