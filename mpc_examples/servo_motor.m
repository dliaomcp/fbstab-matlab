% This file is part of the fbstab-matlab repository
% This project is licensed under the BSD3 license,
% you should have obtained a copy of the license along with this code


clear all;
close all
clc;


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
N = 10;

mpc.N = N;
mpc.nx = nx;
mpc.nu = nu;

%% Cost function *************************************
Q = diag([qy,0,0,0]);
mpc.Q = repmat(Q,[1,1,N+1]);

R = [ru];
mpc.R = repmat(R,[1,1,N+1]);


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

qp = form_mpc_qp_condensed(mpc);

%% Simulate unconstrained control
ytrg = deg2rad(30);

xtrg = [ytrg;0;0;0];
utrg = 0;

x0 = zeros(nx,1);

tfinal = 2;

simin.Tf = tfinal;
simin.ts = ts;
simin.x0 = x0;
simin.xbar = xtrg;
simin.ubar = utrg;

[t,U,X,Y] = simLTI(sys,simin,qp);

figure();
o = ones(size(t));
ax(1) = subplot(3,1,1);
plot(t,rad2deg(Y(1,:)),t,o*rad2deg(ytrg),'--');
legend('signal','target');
ylabel('Shaft Position [\circ]');
ax(2) = subplot(3,1,2);
plot(t,Y(2,:),t,o*ymax,'k--',t,-o*ymax,'k--');
ylabel('Shaft Torque [Nm]');
ax(3) = subplot(3,1,3);
h = stairs(t,U);
hold on 
plot(t,umax*o,'k--',t,-umax*o,'k--');
for i =1:length(h)
	h(i).LineWidth = 1.5;
end
ylabel('Input Voltage [V]')
linkaxes(ax,'x');
xlim([0,tfinal]);
xlabel('Time [s]');












