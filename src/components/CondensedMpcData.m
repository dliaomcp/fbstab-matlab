% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

classdef CondensedMpcData < handle
properties(Access = public)
	N; % prediction horizon
	nx; % state dimension
	nu; % control dimension
	nc; % constraints per stage

	nz; % number of primal variables
	nl; % number of equality constraints
	nv; % number of inequality constraints

	Szero = false; % S = 0
	Hdiag = false; % S = 0 and Q,R are diagonal

	% cost function
	Q; % Q0,Q1,... \in [nx,nx,N+1]
	R; % \in [nu,nu,N+1]
	S; % \in [nu,nx,N+1]
	q; % \in [nx,1]
	r; % \in [nx,1]

	% dynamics
	Ak; % A0,A1,.. AN-1 in [nx,nx,N]
	Bk; % \in [nx,nu,N]
	ck; % in [nx,N]
	xt; % in [nx,1] the state

	% constraints
	E; % \in [nc,nx,N+1]
	L; %\in [nc,nu,N+1]
	d; % \in [nc,N+1]

	% useful vectors
	f;
	b;
end % properties

methods(Access = public)
	% constructor
	function s = CondensedMpcData(Q,R,S,q,r,A,B,c,xt,E,L,d)
		[s.nx,s.nu,s.N] = size(B);
		[nx,nu,N] = size(B);
		s.nc = size(E,1);
		nc = size(E,1);
		% problem sizes
		s.nz = (s.nu)*(s.N+1);
		s.nl = 0;
		s.nv = s.nc*(s.N+1);

		s.Q = Q; s.R = R; s.q = q; s.r = r;
		s.Ak = A; s.Bk = B; s.ck = c; s.E = E;
		s.L = L; s.d = d; s.S = S; s.xt = xt;
		s.Szero = false;
		s.Hdiag = false;

		% form gradient vector
		s.f = zeros((nx+nu)*(N+1),1);
		w = s.evalw();
		r = r(:);
		q = q(:);

		s.f = s.QP(w) + q;
		s.f = s.MT(s.f) + s.SP(w) + r;

		% constraint RHS
		s.b = zeros(nc*(N+1),1);
		s.b = -(d(:) + s.EP(w));

	end

	function [nz,nl,nv] = ProblemSize(obj)
		[nx,nu,nc,N] = OcpSize(obj);
		nz = (nu)*(N+1);
		nl = 0;
		nv = nc*(N+1);
	end

	% compute H*v 
	function y = H(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		% H*v = (R + M'QM + SM + (SM)')*v
		x = obj.MP(v);
		y = obj.RP(v) + obj.SP(x);
		x = obj.QP(x)+obj.ST(v);
		y = y+ obj.MT(x);
		
	end

	% compute A*v
	function y = A(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		% A = EM + L
		y = zeros(nc*(N+1),1);
		y = obj.LP(v);
		y = y + obj.EP(obj.MP(v));
		
	end

	% compute A'*v
	function y = AT(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nu*(N+1));
		y = obj.LT(v) + obj.MT(obj.ET(v));
	
	end

	function [nx,nu,nc,N] = OcpSize(obj)
		[nx,nu,N] = size(obj.Bk);
		nc = size(obj.E,1);
	end

	% compute H*v for the multiple shooting H 
	function y = H_ms(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		% v must be (nx+nu)*(N+1)
		v = reshape(v,[nx+nu,N+1]);
		y = zeros(nx+nu,N+1);
		for i = 1:N+1
			x = v(1:nx,i);
			u = v(nx+1:nx+nu,i);
			y(:,i) = ...
			[obj.Q(:,:,i)*x + obj.S(:,:,i)'*u;
			obj.S(:,:,i)*x + obj.R(:,:,i)*u];
		end
		y = y(:);
	end

	% compute A'*v for multiple shooting A
	function y = AT_ms(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx+nu,N+1);
		v = reshape(v,[nc,N+1]);
		for i = 1:N+1
			y(:,i) = [obj.E(:,:,i)';obj.L(:,:,i)']*v(:,i);
		end
		y = y(:);
	end

	function f = f_ms(obj)
		[nx,nu,nc,N] = OcpSize(obj);
		f = zeros(nx+nu,N+1);

		f(:,1) = [obj.q(:,1);obj.r(:,1)];
		for i = 2:N+1
			f(:,i) = [obj.q(:,i);obj.r(:,i)];
		end
		f = f(:);
	end

end % public methods

methods(Access = {?ImplicitConjugateGradient})
	% Evaluate the w term in x = Mu + w
	function w = evalw(obj)
		[nx,nu,nc,N] = OcpSize(obj);
		w = zeros(nx,N+1);
		w(:,1) = obj.xt;
		for i = 1:N
			w(:,i+1) = obj.Ak(:,:,i)*w(:,i) + obj.ck(:,i);
		end
		w = w(:);
	end

	% State propagation
	function x = MP(obj,u)
		[nx,nu,nc,N] = OcpSize(obj);
		x = zeros(nx,N+1);
		
		uu = reshape(u,[nu,N+1]);
		for i = 1:N
			x(:,i+1) = obj.Ak(:,:,i)*x(:,i) + ...
			obj.Bk(:,:,i)*uu(:,i);
		end
		x = x(:);
	end

	% Co-state propagation
	function y = MT(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nu,N+1);
		v = reshape(v,[nx,N+1]);
		l = zeros(nx,N);
		l(:,N) = v(:,N+1);
		for i = N:-1:2
			l(:,i-1) = v(:,i) + obj.Ak(:,:,i)'*l(:,i);
		end
		for i = 1:N
			y(:,i) = obj.Bk(:,:,i)'*l(:,i);
		end
		y = y(:);
	end

	% Control penalty matrix
	function y = RP(obj,u)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nu,N+1);
		u = reshape(u,[nu,N+1]);
		for i = 1:N+1
			y(:,i) = obj.R(:,:,i)*u(:,i);
		end
		y = y(:);
	end

	% State penalty matrix
	function y = QP(obj,x)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx,N+1);
		x = reshape(x,[nx,N+1]);
		for i = 1:N+1
			y(:,i) = obj.Q(:,:,i)*x(:,i);
		end
		y = y(:);
	end

	% Coupling matrix
	function y = SP(obj,x)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nu,N+1);
		x = reshape(x,[nx,N+1]);
		for i = 1:N+1
			y(:,i) = obj.S(:,:,i)*x(:,i);
		end
		y = y(:);
	end

	function y = ST(obj,u)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx,N+1);
		u = reshape(u,[nu,N+1]);
		for i = 1:N+1
			y(:,i) = obj.S(:,:,i)'*u(:,i);
		end
		y = y(:);
	end

	% State constraint matrix
	function y = EP(obj,x)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nc,N+1);
		x = reshape(x,[nx,N+1]);
		for i = 1:N+1
			y(:,i) = obj.E(:,:,i)*x(:,i);
		end
		y = y(:);
	end

	function y = ET(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx,N+1);
		v = reshape(v,[nc,N+1]);
		for i = 1:N+1
			y(:,i) = obj.E(:,:,i)'*v(:,i);
		end
		y = y(:);
	end

	% Control constraint matrix
	function y = LP(obj,u)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nc,N+1);
		u = reshape(u,[nu,N+1]);
		for i = 1:N+1
			y(:,i) = obj.L(:,:,i)*u(:,i);
		end
		y = y(:);
	end

	function y = LT(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nu,N+1);
		v = reshape(v,[nc,N+1]);
		for i = 1:N+1
			y(:,i) = obj.L(:,:,i)'*v(:,i);
		end
		y = y(:);
	end

end %methods


end