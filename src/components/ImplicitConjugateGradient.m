% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

classdef ImplicitConjugateGradient < handle
	properties(Access = public)
	ztol = 100*eps; % zero tolerance 
	alpha = 0.95;
	pcg_tol = 1e-8; % relative tolerance for PCG iteration
	data;
	sigma;
	mus;
	gamma;
	T; % preconditioner

	end % properties

	methods(Access = public)
	function o = ImplicitConjugateGradient(data)
		o.data = data;
		[nx,nu,nc,N] = data.OcpSize();
		[nz,~,nv] = data.ProblemSize();
		o.gamma = zeros(nv,1);
		o.mus = zeros(nv,1);
		o.T = zeros(nu,N+1);
	end

	function Factor(o,x,xbar,sigma)
		% form the augmented Hessian matrix
		[nx,nu,nc,N] = OcpSize(o.data);
		[nz,nl,nv] = ProblemSize(o.data);

		o.sigma = sigma;
		% compute the barrier terms
		ys = x.y + sigma*(x.v - xbar.v);
		[gamma,mu] = o.dphi(ys,x.v);
		mus = sigma*gamma + mu;
		o.mus = mus;
		o.gamma = gamma;

		for i = 1:N+1
			o.T(:,i) = diag(o.data.R(:,:,i)) + sigma;
		end
	end % factor

	function Solve(o,r,dx)
		[nx,nu,nc,N] = OcpSize(o.data);
		[nz,nl,nv] = ProblemSize(o.data);

		% compute the reduced residual
		r1 = r.rz - o.data.AT(r.rv./o.mus);

		[dx.z,res] = o.ConjugateGradient(r1,1e-8,nz);

		% Recover Inequality duals
		dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus;

		% dy for linesearch
		dx.y = o.data.b - o.data.A(dx.z);
	end

end % public methods

methods(Access = private)
	function [x,res] = ConjugateGradient(o,b,tol,kmax)
		n = length(b);
		x = zeros(n,1);
		r = b;
		bb = dot(b,b);
		tau = dot(r,o.invT(r));
		p = zeros(n,1);
		w = zeros(n,1);
		for k = 1:kmax
			rho = dot(r,r);
			res = sqrt(rho)/bb;
			if res <= tol
				break;
			end
			z = o.invT(r);
			z = r;
			tau1 = tau;
			tau = dot(z,r);
			if k == 1
				p = z;
			else
				beta = tau/tau1;
				p = z+beta*p;
			end
			w = o.Kapply(p);
			alpha = tau/dot(p,w);
			x = x + alpha*p;
			r = r - alpha*w;
		end
	end

	function y = invT(o,u)
		[nx,nu,nc,N] = OcpSize(o.data);
		u = reshape(u,[nu,N+1]);
		y = zeros(nu,N+1);
		for i = 1:N+1
			y(:,i) = u(:,i) ./o.T(:,i);
		end
		y = y(:);
	end

	function y = Kapply(o,u)
		x = o.data.MP(u);
		t = o.data.EP(x) + o.data.LP(u);
		t = t.*o.gamma./o.mus;
		t2 = o.data.QP(x) + o.data.ST(u) + o.data.ET(t);

		y = o.sigma*u + o.data.RP(u) + ...
		o.data.MT(t2) + o.data.SP(x) + o.data.LT(t);
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
	
end % private methods


end % class