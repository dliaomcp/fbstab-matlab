% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 



function [mpc,sys,simin] = load_servo(N)

	if nargin < 1
		N = 10;
	end
	% model parameters
	Ls = 1;
	kt = 10;
	bl = 25;
	Jm = 0.5;
	bm = 0.1;
	ktheta = 1280.2;
	RR = 20;
	rho = 20;
	Jl = 20*Jm;

	% constraints
	umax = 220;
	ymax = 78.5358;

	% sampling
	ts = 0.05;

	% form CT matrices
	Ac = [0,1,0,0;
	-ktheta/Jl,-bl/Jl,ktheta/(rho*Jl),0;
	0,0,0,1;
	ktheta/(rho*Jm),0,-ktheta/(rho^2*Jm),-(bm + kt^2/RR)/Jm];

	Bc = [0;0;0;kt/(RR*Jm)];

	Cc = [1,0,0,0;ktheta,0,-ktheta/rho,0];

	sys = ss(Ac,Bc,Cc,0);
	sysc = sys;
	% convert to discrete time
	sys = c2d(sys,ts,'zoh');
	[nx,nu] = size(sys.B);


	% tuning parameters etc.
	qy = 1000; % tracking penalty
	ru = 1e-4; % control penalty

	%% Cost function *************************************
	Q = diag([qy,0,0,0]);
	mpc.Q = repmat(Q,[1,1,N+1]);

	R = [ru];
	mpc.R = repmat(R,[1,1,N+1]);

	mpc.S = repmat(zeros(nu,nx),[1,1,N+1]);

	%% dynamics *************************************
	mpc.A = repmat(sys.A,[1,1,N]);
	mpc.B = repmat(sys.B,[1,1,N]);
	c = zeros(nx,1);
	mpc.c = repmat(c,[1,N]);


	%% constraints *************************************
	nc = 4;
	mpc.nc = nc;
	% output constraint
	% |y(2)| <= ymax
	E = [sys.C(2,:);-sys.C(2,:);zeros(2,nx)];
	% input constraint
	% |u| <= umax
	L = [0;0;1;-1];

	% can't constraint y(t=0)
	mpc.E = cat(3,zeros(size(E)),repmat(E,[1,1,N]));
	mpc.L = repmat(L,[1,1,N+1]);

	d = [-ymax;-ymax;-umax;-umax];
	mpc.d = repmat(d,[1,N+1]);

	% system object
	ytrg = deg2rad(30);
	xtrg = [ytrg;0;0;0];
	utrg = 0;
	mpc.xtrg = xtrg;
	mpc.utrg = utrg;

	% create q and r terms
	mpc.q = repmat(-Q*xtrg,[1,N+1]);
	mpc.r = repmat(-R*utrg,[1,N+1]);


	sys.name = 'servo';
	simin.x0 = zeros(nx,1);
	simin.tfinal = 2;
	simin.xtrg = xtrg;
	simin.utrg = utrg;
	simin.ytrg = ytrg;
	simin.umax = umax;
	simin.ymax = ymax;
	mpc.x0 = simin.x0;
end