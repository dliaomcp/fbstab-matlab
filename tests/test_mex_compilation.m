% @Author: Dominic Liao-McPherson
% @Date:   2019-03-17 17:07:30
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2019-03-18 00:28:08

% unit tests codegeneration of fbstab_mpc
clear all
close all
clc


[mpc,sys,simin] = load_servo_example(30);

[nx,nu,N] = size(mpc.B);
nc = size(mpc.E,1);

nz = (nx+nu)*(N+1);
nl = nx*(N+1);
nv = nc*(N+1);

% create initial conditions
x.v = zeros(nv,1);
x.l = zeros(nl,1);
x.z = zeros(nz,1);

% options
opts.display_level = 2;
opts.linear_solver = 'ric';

gendir = 'src/bin/';
mkdir(gendir);
addpath(gendir);

args = {x,mpc,opts};
codegen('fbstab_mpc','-config:mex','-o',[gendir 'fbstab_mpc_mex'],'-args',args);


disp('Testing Ricatti liner solver\n\n');
[xopt,out] = fbstab_mpc_mex(x,mpc,opts);

disp('Testing PCG liner solver\n\n');
opts.linear_solver = 'pcg';
[xopt,out] = fbstab_mpc_mex(x,mpc,opts);
