% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


function [mpc,sys,simin] = load_hcw(N)
	if nargin < 1
		N = 40;
	end

	% model parameters
	mu = 398600.441800000; % km^3/s^2
	Re = 6371; % radius of the earth km
	alt = 650; % A reasonable altitude for a LEO orbit km

	% orbital frequency
	n = sqrt(mu/(Re+alt)^3);


	% continuous time model
	A21 = [3*n^2,0,0;0,0,0;0,0,-n^2];
	A22 = [0,2*n,0;-2*n,0,0;0,0,0];
	Ac = [zeros(3),eye(3);A21,A22];
	Bc = [zeros(3);eye(3)];
	C = eye(6);

	sys = ss(Ac,Bc,C,0);

	% transform to discrete
	ts = 30; % sampling period
	sys = c2d(sys,ts,'zoh');

	% form system with delta v as the input
	A = sys.A;
	B = A*[zeros(3);eye(3)];
	sys.A = A;
	sys.B = B;
	C = sys.C;
	[nx,nu] = size(B);

	% MPC params
	qx = 1; % position weight
	qv = 1e-3; % velocity weight
	ru = 1; % control penalty

	umax = 1e-3; %km/s^2
	vmax = 1e-3; % km/s

	%% Cost function *************************************
	Q = diag([qx*ones(3,1);qv*ones(3,1)]);
	mpc.Q = repmat(Q,[1,1,N+1]);

	R = ru*eye(3);
	mpc.R = repmat(R,[1,1,N+1]);

	mpc.S = repmat(zeros(nu,nx),[1,1,N+1]);

	%% Dynamics *************************************
	mpc.A = repmat(A,[1,1,N]);
	mpc.B = repmat(B,[1,1,N]);
	c = zeros(nx,1);
	mpc.c = repmat(c,[1,N]);

	%% Constraints
	% input constraints 
	% |u| <= umax
	E = zeros(2*nu,nx);
	L = [eye(nu);-eye(nu)];
	d = [-umax*ones(3,1);-umax*ones(3,1)];

	% velocity constraints
	% |v| <= vmax
	E = [E;zeros(3),eye(3);-zeros(3),-eye(3)];
	L = [L;zeros(6,nu)];
	d = [d;-vmax*ones(3,1);-vmax*ones(3,1)];

	nc = size(E,1);
	mpc.E = cat(3,zeros(size(E)),repmat(E,[1,1,N]));
	mpc.L = repmat(L,[1,1,N+1]);
	mpc.d = repmat(d,[1,N+1]);


	xtrg = zeros(6,1);
	utrg = zeros(3,1);
	% create q and r terms
	mpc.q = repmat(-Q*xtrg,[1,N+1]);
	mpc.r = repmat(-R*utrg,[1,N+1]);

	mpc.xtrg = xtrg;
	mpc.utrg = utrg;
	
	sys.name = 'hcw';

	%% Simulate
	simin.x0 = [-2.8;-0.01;-1;0;0;0];
	simin.xtrg = zeros(6,1);
	simin.utrg = zeros(3,1);
	simin.vmax = vmax;
	simin.umax = umax;
	simin.tfinal = 3000;
	mpc.x0 = simin.x0;

end

