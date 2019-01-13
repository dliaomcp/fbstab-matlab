% @Author: Dominic Liao-McPherson
% @Date:   2018-10-03 14:54:10
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2018-10-08 16:59:46


% takes in an MPC structure which contains
% N, nx, nu, nc
% A \in [nx,nx,N] = [A0,A1...,A_N-1]
% B \in [nx,nu,N] = [B0,..., B_N-1]
% Q \in [nx,nx,N+1] = [Q0... QN]
% R \in [nx,nx,N+1]
% E \in [nc,nx,N+1]
% L \in [nc,nu,N+1]
% xbar \in [nx,N+1]
% ubar \in [nu,N+1]
% c \in [nx,N]
% d \in [nx,N+1]
% 
% Which define the QP
% min   \sum  ||x(i)-xbar||_Q^2 + ||u(i)-ubar||_R^2
% s.t., x(i+1) =  A(i)x(i)+B(i)u(i) +c(i),  x(0) = x(t)
%       E(i)x(i) + L(i)u(i) + d(i) <= 0
%
% x = [x0;x1;...xN],  u = [u0;u1;...uN]
%
% outputs the matrices necessary to form
% min 1/2 u'Hu + (h + Tx xbar + Tu ubar + Tf x(t))
% s.t.  G u <=  g + Tb x(t)
%

% outputs the following 
function qp = form_mpc_qp(mpc)
	NN = mpc.N+1;
	nx = mpc.nx;
	nu = mpc.nu;
	nc = mpc.nc;
	N = mpc.N;

	% get the state transition matrices
	[qp.M,qp.T,qp.w] = stm(mpc);
	M = qp.M;

	%% Cost function formation
	% form the big Q and R matrices
	R = zeros(nu*NN);
	Q = zeros(nx*NN);
	for i = 1:NN
		im = 1 + (i-1)*nx;
		ip = i*nx;
		Q(im:ip,im:ip) = mpc.Q(:,:,i);

		im = 1+(i-1)*nu;
		ip = i*nu;
		R(im:ip,im:ip) = mpc.R(:,:,i);
	end
	qp.R = R;
	qp.Q = Q;

	% measurement input matrix
	qp.Tf = M'*Q*qp.T;

	% reference input matrix
	II = kron(ones(NN,1),eye(nx));
	qp.Tx = -M'*Q*II;

	% nominal control input matrix
	II = kron(ones(NN,1),eye(nu));
	qp.Tu = -R*II;

	% constant term
	qp.h = M'*Q*qp.w;

	% reduced Hessian
	qp.H = R + M'*Q*M;

	%% Constraint formation
	% form the big E and L matrices
	E = zeros(NN*nc,NN*nx);
	L = zeros(NN*nc,NN*nu);
	for i = 1:NN
		im = 1 + nc*(i-1);
		ip = nc*i;
		jm = 1+nx*(i-1);
		jp = nx*i;

		E(im:ip,jm:jp) = mpc.E(:,:,i);

		jm = 1+nu*(i-1);
		jp = nu*i;

		L(im:ip,jm:jp) = mpc.L(:,:,i);
	end
	qp.E = E;
	qp.L = L;

	% vector portion
	d = mpc.d(:);
	qp.g = -(d + E*qp.w);
	% measurement input matrix
	qp.Tb = -E*qp.T;
	% constraint Jacobian
	qp.G = E*M+L;

	% sizes

	qp.nx = nx;
	qp.nu = nu;
	qp.nz = size(qp.H,1);
	qp.N = N;
	qp.nv = size(qp.G,1);
	qp.nl = 0;
end


% return the matrices such that
% x = Mu + T x(t) + w
function [M,T,w] = stm(s)
	n = s.nx;
	m = s.nu;
	N = s.N;
	M = zeros(n*(N+1),m*(N+1));
	w = zeros(n,(N+1));

	Ap = zeros(n,n,N);
	Ap(:,:,1) = eye(n);
	for i = 2:N
		Ap(:,:,i) = s.A(:,:,i)*Ap(:,:,i-1);
	end

	% form the matrix
	for i = 1:N
		for k = 1:N
			if i >= k
				im = n*(i-1)+1 + n;
				ip = i*n + n;
				km = m*(k-1)+1;
				kp = k*m;
				idx = max(1,i-k+1);
				M(im:ip,km:kp) = Ap(:,:,idx)*s.B(:,:,k);
			end
		end
	end
	% form the IC input gain
	% w = [w0,w1...wN]
	% start with w0 = 0, w1 = c0
	w(:,2) = s.c(:,1);
	for k = 3:N+1
		% w(k) = A(k-1)*w(k-1) + c(k)
		w(:,k) = s.A(:,:,k-1)*w(:,k-1) + s.c(:,k-1);
	end
	w = w(:);

	% form the MPC state gain
	% T = [I;A0;A1*A0;A2*A1*A0...]
	T = cat(3,eye(n),repmat(s.A(:,:,1),[1,1,N]));
	for k = 2:N+1
		T(:,:,k) = Ap(:,:,k-1)*T(:,:,k);
	end
	T = permute(T,[1,3,2]);
	T = reshape(T,[(N+1)*n,n]);
end