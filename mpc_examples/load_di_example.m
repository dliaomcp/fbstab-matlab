% @Author: Dominic Liao-McPherson
% @Date:   2019-03-13 16:15:54
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2019-03-13 16:16:22


function mpc = load_di(N)

	% simple double integrator example
	Q = diag([2,1]);
	S = [1,0];
	R = 3;
	xtrg = [1;0];
	q = -Q*xtrg;
	r = 0;

	A = [1,1;0,1];
	B = [0;1];
	c = [0;0];
	x0 = [0;0];

	E = [-eye(2);eye(2);zeros(2,2)];
	L = [zeros(4,1);-1;1];
	d = [0;0;-2;-2;-1;-1];

	% create the data structure
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