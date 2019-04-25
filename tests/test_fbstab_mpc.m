% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

clear all
close all
clc

% call the servo motor example
[mpc,sys,simin] = load_servo_example(50);

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

'MA57'
tic
[xopt,out] = fbstab_ldl_sp(x,mpc,opts);
toc



