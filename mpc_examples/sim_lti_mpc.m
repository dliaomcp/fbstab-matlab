% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 



function [T,U,X,Y] = sim_lti_mpc(sys,simin,mpc)

x0 = simin.x0;
tf = simin.tfinal;
xtrg = simin.xtrg;
utrg = simin.utrg;

ts = sys.ts;
M = ceil(tf/ts);

A = sys.A;
B = sys.B;
C = sys.C;

[nx,nu,N] = size(mpc.B);
nc = size(mpc.E,1);
ny = size(C,1);

X = zeros(nx,N+1);
U = zeros(nu,N+1);
Y = zeros(ny,N+1);
T = zeros(N+1,1);

X(:,1) = x0;
Y(:,1) = C*x0;

s.z = zeros((nx+nu)*(N+1),1);
s.l = zeros(nx*(N+1),1);
s.v = zeros(nc*(N+1),1);

t = 0;
for i = 1:M
	x = X(:,i);
	[u,s] = lmpc(s,x,mpc,xtrg,utrg);
	U(:,i) = u;

	X(:,i+1) = A*x + B*u;
	Y(:,i+1) = C*X(:,i+1);

	T(i) = t;
	t = t+ts;
end
T(end+1) = t;
U(:,end+1) = U(:,end);


end

function [u,s] = lmpc(s,x,mpc,xtrg,utrg)
	[nx,nu,N] = size(mpc.B);

	mpc.x0 = x;
	for i = 1:N+1
		mpc.q(:,i) = -mpc.Q(:,:,i)*xtrg - mpc.S(:,:,i)'*utrg;
		mpc.r(:,i) = -mpc.R(:,:,i)*utrg - mpc.S(:,:,i)*xtrg;
	end

	[s,u] = fbstab_mpc(s,mpc);
end

