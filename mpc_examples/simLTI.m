% @Author: Dominic Liao-McPherson
% @Date:   2018-10-08 15:28:37
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2019-01-13 15:10:25


%% performs a closed-loop simulation of an LTI model
% u = ctrl(t,x)
function [T,U,X,Y] = simLTI(sys,simin,qp)

	A = sys.A;
	B = sys.B;
	C = sys.C;

	[nx,nu] = size(B);
	ny = size(C,1);

	N = ceil(simin.Tf/simin.ts);
	X = zeros(nx,N+1);
	U = zeros(nu,N+1);
	Y = zeros(ny,N+1);
	T = zeros(N+1,1);
	X(:,1) = simin.x0;
	Y(:,1) = C*simin.x0;

	s.z = zeros(qp.nz,1);
	s.v = zeros(qp.nv,1);

	t = 0;
	for i = 1:N
		x = X(:,i);
		[u,s,res] = lmpc(x,simin.xbar,simin.ubar,qp,s);

		U(:,i) = u;

		X(:,i+1) = A*x + B*u;
		Y(:,i+1) = C*X(:,i+1);

		T(i) = t;
		t = t+simin.ts;
	end
	T(end) = t;
	U(:,end) = U(:,end-1);
end