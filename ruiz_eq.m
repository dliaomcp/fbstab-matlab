% @Author: dliaomcp
% @Date:   2018-04-25 19:02:21
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2018-08-14 14:19:58
% 
% [sc,qp] = ruiz_equlibration(qp,opts)
% Compute a diagonal scaling and cost scaling c which 
% reduces the condition number of M = [H,A';A,0]
% The original problem is 
% min. 1/2 x'Hx + f'x
% s.t. Ax <= b
% initial conditions should be rescaled as
% x0 = x0./p, v0 = v0./d * c
% the solutions can be recovered as
% xopt = xopt.*p, vopt = 1/c * vopt.*d


function [sc,qp,out] = ruiz_eq(qp,opts)
	if nargin < 2
		opts = struct();
	end

	max_ruiz_iters = 30;
	ruiz_tol = 1e-2;
	scale_cost = false;

	if isfield(opts,'max_ruiz_iters')
		max_ruiz_iters = opts.max_ruiz_iters;
	end
	if isfield(opts,'ruiz_tol')
		ruiz_tol = opts.ruiz_tol;
	end
	if isfield(opts,'scale_cost')
		scale_cost = opts.scale_cost;
	end

	n = length(qp.f);
	q = length(qp.b);
	N = n+q;

	c = 1;
	delta = zeros(N,1);
	s = ones(N,1);
	gamma = 1;

	M = [qp.H,qp.A';qp.A,zeros(q)];
	temp = zeros(n,1);
	j = 0;
	while infnorm(1-delta) > ruiz_tol && j <= max_ruiz_iters
		for i = 1:N
			rownorm = infnorm(M(i,:));
			% safeguard against large scalings
			if rownorm < 1e-4
				rownorm = 1;
			end
			delta(i) = 1/sqrt(rownorm);
		end
		M = lrscale(M,delta,delta);
		qp.f = delta(1:n).*qp.f;
		qp.b = delta(n+1:end).*qp.b;

		if scale_cost
			for i = 1:n
				temp(i) = infnorm(M(i,1:n));
			end
			cnorm = max(mean(temp),infnorm(qp.f));
			% safeguard
			if cnorm < 1e-4
				cnorm = 1;
			end
			gamma = 1/cnorm;
			M(1:n,1:n) = gamma*M(1:n,1:n);
			qp.f = gamma*qp.f;
			c = gamma*c;
		end
		s = delta.*s;
		j = j+1;
	end
	qp.H = M(1:n,1:n);
	qp.A = M(n+1:end,1:n);
	p = s(1:n);
	d = s(n+1:end);

	sc.pr = p; % primal scaling
	sc.ieq = d; % dual scaling
	sc.c = c; % cost scaling
	out.iters = j;
end

function y = infnorm(v)
	y = max(abs(v));
end

% computes A = diag(l)*A*diag(r)
function A = lrscale(A,l,r)
	[m,n] = size(A);

	% scale the rows by l
	for i = 1:m
		A(i,:) = l(i)*A(i,:);
	end

	% scale the columns by r
	for i = 1:n
		A(:,i) = r(i)*A(:,i);
	end
end

