% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% Transforms the problem data for the following
% optimal control problem
%
% min.  \sum_{i=0}^N 1/2 [x(i)]' * [Q(i) S(i)'] [x(i)] + [q(i)]'*[x(i)]
%                        [u(i)]    [S(i) R(i) ] [u(i)]   [r(i)]  [u(i)]
% s.t.  x(i+1) = A(i)*x(i) + B(i) u(i) + c(i), i = 0 ... N-1
%       x(0) = x0
%       E(i)*x(i) + L(i)*u(i) + d(i) <= 0,     i = 0 ... N
%
% where [Q(i),S(i)';S(i),R(i)] is positive semidefinite, and
% transforms it into the following form:
%
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az ≤ b,
% 
% with all the data saved as sparse matrices.
%

function qp = form_sparse_mpc_qp(mpc)
	[nx,nu,N] = size(mpc.B);
	nc = size(mpc.E,1);
	NN = N+1;

	nz = NN*(nx+nu);
	nl = NN*nx;
	nv = nc*NN;

	%% Cost function formation
	% form the big Q and R matrices
	qp.H = sparse(nz,nz);

	for i = 1:N+1
		im = (i-1)*(nx+nu) + 1;
		ip = i*(nx+nu);
		qp.H(im:ip,im:ip) = ...
		[mpc.Q(:,:,i),mpc.S(:,:,i)';mpc.S(:,:,i),mpc.R(:,:,i)];
	end

	% G matrix
	qp.G = sparse(nl,nz);
	qp.G(1:nx,1:nx) = -eye(nx);
	for i = 2:N+1
		im = (i-1)*nx+1;
		ip = i*nx;

		jm = (i-2)*(nx+nu)+1;
		jp = jm + 2*nx+nu-1;

		qp.G(im:ip,jm:jp) = ...
		[mpc.A(:,:,i-1),mpc.B(:,:,i-1),-eye(nx)];
	end

	% A matrix
	qp.A = sparse(nv,nz);
	for i = 1:N+1
		im = 1+(i-1)*nc;
		ip = i*nc;

		jm = 1+(i-1)*(nx+nu);
		jp = i*(nx+nu);

		qp.A(im:ip,jm:jp) = [mpc.E(:,:,i),mpc.L(:,:,i)];
	end


	qp.f = reshape([mpc.q;mpc.r],[nz,1]);
	qp.h = -[mpc.x0;mpc.c(:)];
	qp.b = -mpc.d(:);
end
