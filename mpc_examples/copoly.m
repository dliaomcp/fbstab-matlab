% @Author: Dominic Liao-McPherson
% @Date:   2018-09-21 12:22:37
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2019-01-13 15:11:34


clear all
close all
clc

% construct transfer matrix

s = tf('s');

G(1,1) = 0.34/(0.85*s+1);
G(1,2) = 0.21/(0.42*s+1);
G(1,3) = 0.5*(0.5*s+1)/(12*s^2+0.4*s+1);
G(1,5) = 6.46*(0.9*s+1)/(0.007*s^2 + 0.3*s + 1);

G(2,1) = -0.41/(2.41*s+1);
G(2,2) = 0.66/(1.51*s+1);
G(2,3) = -0.3/(1.45*s+1);
G(2,5) = -0.372/(0.8*s+1);

G(3,1) = 0.3/(2.54*s+1);
G(3,2) = 0.49/(1.54*s+1);
G(3,3) = -0.71/(1.35*s+1);
G(3,4) = -0.2/(2.71*s+1);
G(3,5) = -4.71/(0.008*s^2 + 0.41*s+1);
G(4,5) = 1.02*(0.23*s+1)/(0.07*s^2 + 0.31*s +1);


% Form a realization of the system
sys = ss(G);

% convert to discrete time
ts = 0.5;
sys = c2d(sys,ts,'zoh');

[nx,nu] = size(sys.B);
ny = size(sys.C,1);
% tuning parameters
N = 70;
qy = 1;
ru = 0.1;
mpc.N = N;
mpc.nx = nx;
mpc.nu = nu;

%% cost function*************************************
Q = qy*(sys.C)'*sys.C;
R = ru*eye(nu);
Qf = dlyap(sys.A,Q);

mpc.Q = repmat(Q,[1,1,N]);
mpc.Q = cat(3,mpc.Q,Qf);
mpc.R = repmat(R,[1,1,N+1]);

%% dynamics *************************************
mpc.A = repmat(sys.A,[1,1,N]);
mpc.B = repmat(sys.B,[1,1,N]);
c = zeros(nx,1);
mpc.c = repmat(c,[1,N]);

%% Constraints
% |u| <= umax
umax = 5/100;
nc = 2*nu;
E = zeros(nc,nx);
L = [eye(nu);-eye(nu)];

mpc.E = repmat(E,[1,1,N+1]);
mpc.L = repmat(L,[1,1,N+1]);
d = [-ones(nu,1)*umax;-ones(nu,1)*umax];
mpc.d = repmat(d,[1,N+1]);
mpc.nc = nc;

% form QP matrices

qp = form_mpc_qp_condensed(mpc);

%% Simulation
x0 = zeros(nx,1);
for j = 1:nx
	x0(j) = 0.2*sin(j);
end

xtrg = zeros(nx,1);
utrg = zeros(nu,1);

tfinal = 100;

simin.Tf = tfinal;
simin.ts = ts;
simin.x0 = x0;
simin.xbar = xtrg;
simin.ubar = utrg;

[t,U,X,Y] = simLTI(sys,simin,qp);

figure();
subplot(2,1,1)
plot(3*t,Y);
legend('y_1','y_2','y_3','y_4','orientation','horizontal');
ylabel('Normalized Outputs');

o = ones(size(t));
subplot(2,1,2);
plot(3*t,U,3*t,umax*o,'k--',3*t,-umax*o,'k--');
legend('u_1','u_2','u_3','u_4','u_5','orientation','horizontal');
ylabel('Normalized Inputs');
xlabel('Time [hr]');






