% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

clear all
close all
clc


Q = diag([0;0]);
R = 0;
S = [0,0];
q = [-1;-1];
r = 0;

% Double integrator matrices
% x = [r;v]
A = [1,1;0,1];
B = [0;1];
c = [0;0];
x0 = [0;0];

nx = 2;
nu = 1;

fprintf('Degenerate Double Integrator Example:\n\n');
N = 4;

E = [1,0;-1,0;0,0;0,0;0,0];
L = [0;0;1;-1;0];
d = [-4;4;-1;-1;0];
nc = length(d);

mpc = ExtendOverHorizon(N,Q,R,S,q,r,A,B,c,x0,E,L,d);
for i = 1:N
  mpc.E(:,:,i) = zeros(nc,nx);
  mpc.d(1:2,i) = [0;0];
end

nz = (nx+nu)*(N+1);
nl = nx*(N+1);
nv = nc*(N+1);
x.z = zeros(nz,1);
x.l = zeros(nl,1);
x.v = zeros(nv,1);

[x,~,out] = fbstab_mpc(x,mpc);

fprintf('Infeasible Double Integrator Example:\n\n');
N = 3;

% r(N) = 4 -> 4 <= r(N) <= 4
% |u(i)| <= 1
E = [1,0;-1,0;0,0;0,0];
L = [0;0;1;-1];
d = [-4;4;-1;-1];
nc = length(d);

mpc = ExtendOverHorizon(N,Q,R,S,q,r,A,B,c,x0,E,L,d);
for i = 1:N
  mpc.E(:,:,i) = zeros(nc,nx);
  mpc.d(1:2,i) = [0;0];
end

nz = (nx+nu)*(N+1);
nl = nx*(N+1);
nv = nc*(N+1);
x.z = zeros(nz,1);
x.l = zeros(nl,1);
x.v = zeros(nv,1);

[x,~,out] = fbstab_mpc(x,mpc);

fprintf('Unbounded Below Double Integrator Example:\n\n');
N = 3;
Q = diag([0;0]);
R = 0;
S = [0,0];
q = [-1;-1];
r = 0;

% u(i) >= 0 
E = [0,0];
L = [-1];
d = [0];
nc = length(d);
mpc = ExtendOverHorizon(N,Q,R,S,q,r,A,B,c,x0,E,L,d);

nz = (nx+nu)*(N+1);
nl = nx*(N+1);
nv = nc*(N+1);
x.z = zeros(nz,1);
x.l = zeros(nl,1);
x.v = zeros(nv,1);

[x,~,out] = fbstab_mpc(x,mpc);


function mpc = ExtendOverHorizon(N,Q,R,S,q,r,A,B,c,x0,E,L,d)
  mpc.Q = repmat(Q,[1,1,N+1]);
  mpc.R = repmat(R,[1,1,N+1]);
  mpc.S = repmat(S,[1,1,N+1]);
  mpc.q = repmat(q,[1,N+1]);
  mpc.r = repmat(r,[1,N+1]);

  mpc.A = repmat(A,[1,1,N]);
  mpc.B = repmat(B,[1,1,N]);
  mpc.c = repmat(c,[1,N]);
  mpc.x0 = x0;

  mpc.E = repmat(E,[1,1,N+1]);
  mpc.L = repmat(L,[1,1,N+1]);
  mpc.d = repmat(d,[1,N+1]);
end