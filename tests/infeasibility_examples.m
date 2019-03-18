% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

clear all
close all
clc

% script to test the fbstab solver on degenerate and primal/dual infeasible QPs

x0 = [0;0];
opts = struct();

% degenerate
fprintf('Degenerate but feasible QP, eflag should be 0 and dx, dv should be near zero\n');
H = [1,0;0,0];
f = [1;0];
A = [0,0;1,0;0,1;-1,0;0,-1];
b = [0;3;3;-1;-1];
qp.H = H;
qp.f = f;
qp.A = A;
qp.b = b;
v0 = zeros(size(b));
[x,v,out] = fbstab_dense(qp,x0,v0,opts);
out

% infeasible
fprintf('QP is  infeasible, eflag should be -2 and dv should be nonzero\n');
H = [1,0;0,0];
f = [1;-1];
A = [1,1;1,0;0,1;-1,0;0,-1];
b = [0;3;3;-1;-1];
qp.H = H;
qp.f = f;
qp.A = A;
qp.b = b;

v0 = zeros(size(b));
[x,v,out] = fbstab_dense(qp,x0,v0,opts);
out

% unbounded below, [0;1] is an unbounded direction
fprintf('QP is unbounded below, eflag should be -3 and dx should be a unbounded descent direction\n');
H = [1,0;0,0];
f = [1;-1];
A = [0,0;1,0;-1,0;0,-1];
b = [0;3;-1;-1];
qp.H = H;
qp.f = f;
qp.A = A;
qp.b = b;

v0 = zeros(size(b));

[x,v,out] = fbstab_dense(qp,x0,v0,opts);
out





