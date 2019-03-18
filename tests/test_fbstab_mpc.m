% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

clear all
close all
clc

% call the servo motor example
[mpc,qpc,qp,sys,simin] = load_servo(50);

nx = qpc.nx;
nu = qpc.nu;
N = qpc.N;

alpha = 0.95;
% create the data and linsys objects
data = data_ms(mpc.Q,mpc.R,mpc.S,mpc.q,mpc.r,mpc.A,mpc.B,mpc.c,simin.x0,mpc.E,mpc.L,mpc.d);
mpc.x0 = simin.x0;
[nz,nl,nv] = data.opt_sz();

opts.display_level = 2;

% create initial conditions
x.v = zeros(nv,1);
x.l = zeros(nl,1);
x.z = zeros(nz,1);


'implicit pcg'
opts.linear_solver = 'pcg';
tic
[xopt,out] = fbstab_mpc(x,mpc,opts);
toc

'ricatti'
opts.linear_solver = 'ricatti';
tic
[xopt,out] = fbstab_mpc(x,mpc,opts);
toc

% get control sequence
z = xopt.z;
z = reshape(z,[nx+nu,N+1]);
uopt = z(nx+1:end,:);
uopt = uopt(:);
uopt = reshape(uopt,[nu,N+1]);

% compare with the condensed solution
s.H = (qpc.H + qpc.H')/2;
s.f = qpc.h + qpc.Tx*simin.xtrg + qpc.Tu*simin.utrg + qpc.Tf*simin.x0;
s.A = qpc.G;
s.b = qpc.g + qpc.Tb*simin.x0;

'dense quadprog'
tic
u = quadprog(s.H,s.f,s.A,s.b);
toc
u = reshape(u,[nu,N+1]);

should_be_zero = norm(u - uopt)


