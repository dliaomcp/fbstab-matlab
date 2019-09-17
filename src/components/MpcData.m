% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% a class to store the MPC problem data
% this class assumes that the primal variable is ordered as
% x0,u0,x1,u1,... xN,uN
classdef MpcData < handle

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

	% Cost function
	Q; % Q0,Q1,... \in [nx,nx,N+1]
	R; % \in [nu,nu,N+1]
	S; % \in [nu,nx,N+1]
	q; % \in [nx,1]
	r; % \in [nx,1]

	% Dynamics
	Ak; % A0,A1,.. AN-1 in [nx,nx,N]
	Bk; % \in [nx,nu,N]
	ck; % in [nx,N]
	xt; % in [nx,1] the state

	% Constraints
	E; % \in [nc,nx,N+1]
	L; %\in [nc,nu,N+1]
	d; % \in [nc,N+1]

	% Useful vectors
	f;
	h;
	b;

end % properties

methods(Access = public)
	% constructor
	function s = MpcData(Q,R,S,q,r,A,B,c,xt,E,L,d)
		[s.nx,s.nu,s.N] = size(B);
		[nx,nu,N] = size(B);
		s.nc = size(E,1);
		nc = size(E,1);
		% problem sizes
		s.nz = (s.nx+s.nu)*(s.N+1);
		s.nl = s.nx*(s.N+1);
		s.nv = s.nc*(s.N+1);

		s.Q = Q; s.R = R; s.q = q; s.r = r;
		s.Ak = A; s.Bk = B; s.ck = c; s.E = E;
		s.L = L; s.d = d; s.S = S; s.xt = xt;
		s.Szero = false;
		s.Hdiag = false;

		% form the useful vectors
		f = zeros(nx+nu,N+1);
		h = zeros(nx,N+1);
		b = zeros(nc,N+1);

		f(:,1) = [s.q(:,1);s.r(:,1)];
		h(:,1) = -s.xt;
		for i = 2:N+1
			f(:,i) = [s.q(:,i);s.r(:,i)];
			h(:,i) = -s.ck(:,i-1);
		end
		s.f = f(:);
		s.h = h(:);
		s.b = -d(:);
	end

	function [nz,nl,nv] = ProblemSize(obj)
		% The sizes are recomputed to help
		% with code generation.
		[nx,nu,nc,N] = OcpSize(obj);
		nz = (nx+nu)*(N+1);
		nl = nx*(N+1);
		nv = nc*(N+1);
	end

	% compute H*v 
	function y = H(obj,v)
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

	% compute A*v
	function y = A(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nc,N+1);
		v = reshape(v,[nx+nu,N+1]);
		for i = 1:N+1
			x = v(1:nx,i);
			u = v(nx+1:nx+nu,i);
			y(:,i) = obj.E(:,:,i)*x + obj.L(:,:,i)*u;
		end
		y = y(:);
	end

	% compute G*v 
	function y = G(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx,N+1);
		v = reshape(v,[nx+nu,N+1]);
		y(:,1) = -v(1:nx,1); % -x0
		for i = 2:N+1
			xm1 = v(1:nx,i-1);
			um1 = v(nx+1:nx+nu,i-1);
			x = v(1:nx,i);
			y(:,i) = ...
			obj.Ak(:,:,i-1)*xm1 + obj.Bk(:,:,i-1)*um1 - x;
		end
		y = y(:);
	end

	% compute G'*v
	function y = GT(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx+nu,N+1);
		v = reshape(v,[nx,N+1]);
		for i = 1:N
			lm1 = v(:,i);
			l = v(:,i+1);
			y(:,i) = [-lm1 + obj.Ak(:,:,i)'*l;obj.Bk(:,:,i)'*l];
		end
		lN = v(:,N+1);
		y(:,N+1) = [-lN;zeros(nu,1)];
		y = y(:);
	end

	% compute A'*v
	function y = AT(obj,v)
		[nx,nu,nc,N] = OcpSize(obj);
		y = zeros(nx+nu,N+1);
		v = reshape(v,[nc,N+1]);
		for i = 1:N+1
			y(:,i) = [obj.E(:,:,i)';obj.L(:,:,i)']*v(:,i);
		end
		y = y(:);
	end
	
	function [nx,nu,nc,N] = OcpSize(obj)
		% The sizes are recomputed to help
		% with code generation.
		[nx,nu,N] = size(obj.Bk);
		nc = size(obj.E,1);
	end

end %methods


end