% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

classdef fbstab_algo < handle

properties(Access = public)
	% algorithm parameters
	sigma = sqrt(eps);
	max_newton_iters = 500;
	max_prox_iters = 100;
	tol = 1e-6;
	rtol = 1e-12;
	dtol = 1e-12;
	inf_tol = 1e-8;
	check_infeasibility = true;
	alpha = 0.95;
	beta = 0.7;
	eta = 1e-8;
	lsmax = 20;
	max_inner_iters = 100;
	itol_max = 1e-1;
	itol_min = 10*eps;
	itol_red_factor = 1/5;

	% display settings
	% 0 = none
	% 1 = final
	% 2 = outer iters
	% 3 = inner iters
	display_level = 0;

	% constituent classes
	data; % qp data object
	linsys; % linear system object
	feas; % feasibility checker
	% residual objects
	rk; 
	ri;

	% variable objects
	xi;
	xk;
	xp;
	dx;


	% problem size data
	nz;
	nv;
	nl;

end

methods(Access = public)
	% constructor 
	function o = fbstab_algo(data,linsys,feas,rk,ri,xk,xi,xp,dx,opts)
		if nargin < 5
			opts = struct();
		else
			o.update_options(opts);
		end

		o.data = data;
		o.linsys = linsys;
		o.feas = feas;
		o.rk = rk;
		o.ri = ri;

		o.xk = xk;
		o.xi = xi;
		o.xp = xp;
		o.dx = dx;

		% harmonize alpha
		o.rk.alpha = o.alpha;
		o.ri.alpha = o.alpha;
		linsys.alpha = o.alpha;

		[nz,nl,nv] = data.opt_sz();
	end

	% solve the qp
	function [x,out] = solve(o,x0)
		% get problem sizes
		[nz,nl,nv] = opt_sz(o.data);

		%% initialization *************************************
		xk = o.xk; % initial point
		xi = o.xi; % variable for the pfb solver
		dx = o.dx; % change in x and newton step
		xp = o.xp; % linesearch query point

		xk.scopy(x0);
		dx.ones();

		% residual objects
		rk = o.rk; % outer loop
		ri = o.ri; % inner loop
		rk.nres(xk);
		
		% *************************************

		% quality of initial guess
		E0 = rk.norm();
		dE = 1e6;

		% initialize output struct
		out = struct();
		out.eflag = -1;
		out.newton_iters = 0;
		out.prox_iters = 0;
		out.res = E0;

		% initialize tolerance
		itol = min(E0,o.itol_max);
		itol = max(itol,o.itol_min);

		% initial regularization parameter
		sigma = o.sigma;
		
		o.print_iter_header();

		% main proximal loop
		for k = 1:o.max_prox_iters
			% termination check
			rk.nres(xk);
			Ek = rk.norm();

			%% TODO, issue is that sigma parameter needs to vary
			% a large sigma value works for z and v
			% a small value is needed to force a reduction in l
			% I need some kind of adaptive heuristic
			if Ek <= o.tol + E0*o.rtol
				out.eflag = 0;
				break;
			elseif out.newton_iters > o.max_newton_iters
				out.eflag = -1;
				break;
			elseif dE <= o.dtol
				out.eflag = -2; %
			end

			o.print_detailed_header(out.prox_iters,out.newton_iters,rk);
			o.print_iter_line(out.prox_iters,out.newton_iters,rk,ri,itol);
				
			% update inner tolerance
			itol = min(Ek,itol*o.itol_red_factor);
			itol = max(itol,o.itol_min);

			% Call inner solver *************************************
			[Ei,out,flag] = o.eval_prox_pfb(sigma,itol,Ek,out);

			if ~flag
				break;
			end

			% dx = xi - xk	
			dx.vcopy(xi);
			dx.axpy(xk,-1);
			xk.vcopy(xi);
			% infeasibility check
			if o.check_infeasibility
				[primal,dual] = o.feas.check_feas(dx,o.inf_tol);
				if ~primal && ~dual
					out.eflag = -5;
					x = dx.swrite();
					out.res = norm(rk);
					o.print_final(out,rk);
					return;
				elseif ~primal
					out.eflag = -3;
					x = dx.swrite();
					out.res = norm(rk);
					o.print_final(out,rk);
					return;

				elseif ~dual
					out.eflag = -4;
					x = dx.swrite();
					out.res = norm(rk);
					o.print_final(out,rk);
					return;
				end
			end
			out.prox_iters = out.prox_iters +1;
		end % prox loop

		% output solution
		rk.nres(xk);
		out.res = norm(rk);
		x = xk.swrite();
		o.print_final(out,rk);

	end % end solve

	function [Ei,out,flag] = eval_prox_pfb(o,sigma,itol,Epen,out)
		x = o.xi;
		xbar = o.xk;
		dx = o.dx;
		xp = o.xp;
		ri = o.ri;
		rk = o.rk;

		flag = true;
		mrec = zeros(5,1);
		t = 1;

		Ei = 0;
		x.vcopy(xbar);
		for j = 1:o.max_inner_iters
			% compute residuals
			ri.calcres(x,xbar,sigma);
			rk.nres(x);
			% convergence check
			Ei = ri.norm();
			dx.vcopy(x);
			dx.axpy(xbar,-1);
			if (Ei <= itol*min(1,norm(dx)) && norm(rk) < Epen)||(Ei <= o.itol_min) 
				o.print_detailed_footer(itol,ri);
				return
			end

			o.print_detailed_line(j,t,ri);

			if out.newton_iters >= o.max_newton_iters
				out.eflag = -1;
				rk.nres(xbar);
				out.res = norm(rk);
				x = xbar.swrite();
				o.print_final(out,rk);
				flag = false;
				return
			end

			% form and factor the iteration matrix
			o.linsys.factor(x,xbar,sigma);

			% solve the system for the rhs -ri
			% store the result in dx
			ri.negate();
			o.linsys.solve(ri,dx);

			% linesearch
			mrec = circshift(mrec,1);
			mrec(1) = 1/2*Ei^2;
			m0 = max(mrec);
			t = 1;
			for i = 1:o.lsmax
				% compute xp = x + t*dx
				xp.vcopy(x);
				xp.axpy(dx,t);
				ri.calcres(xp,xbar,sigma);
				mp = 1/2*norm(ri)^2;

				if mp <= m0 - 2*t*o.eta*mrec(1)
					break;
				else
					t = o.beta*t;
				end
			end
			% update x = x + t*dx
			x.axpy(dx,t);
			out.newton_iters = out.newton_iters+1;
		end % pfb loop
		% project duals onto the nonnegative orthant
		x.dual_proj();
	end

	% printing
	function print_final(o,out,r)
		if o.display_level >= 1
			switch out.eflag
				case 0
					fprintf('Success\n');
				case -1
					fprintf('Iteration limit exceeded\n');
				case -2
					fprintf('Algorithm stalled\n');
				case -3
					fprintf('Infeasibility detected\n');
				case -4
					fprintf('Unboundedness detected\n');
				otherwise
					fprintf('????');
			end

			[rz,rl,rv] = r.norms();

			a1 = int32(out.prox_iters);
			a2 = int32(o.max_prox_iters);

			fprintf('Proximal iterations: %d out of %d\n', a1, a2);

			a3 = int32(out.newton_iters);
			a4 = int32(o.max_newton_iters);
			fprintf('Newton iterations: %d out of %d\n', a3, a4);
			fprintf('%10s  %10s  %10s  %10s\n','|rz|','|rl|','|rv|','Tolerance');
			fprintf('%10.4e  %10.4e  %10.4e  %10.4e\n',rz,rl,rv,o.tol);
		end
	end

	function print_iter_header(o)
		if o.display_level == 2
			fprintf('%12s %12s %12s %12s %12s %12s %12s \n','prox iter','newton iters','|rz|','|rl|','|rv|','Inner res','Inner tol');
		end
	end

	function print_iter_line(o,prox_iters,newton_iters,rk,ri,itol)
		if o.display_level == 2
			[rz,rl,rv] = rk.norms();
			a1 = int32(prox_iters);
			a2 = int32(newton_iters);
			fprintf('%12d %12d %12.4e %12.4e %12.4e %12.4e %12.4e\n',a1,a2,rz,rl,rv,ri.norm(),itol);
		end
	end

	function print_detailed_header(o,prox_iters, newton_iters, r)
		if o.display_level == 3
			a1 = int32(prox_iters);
			a2 = int32(newton_iters);
			fprintf('Begin Prox Iter: %d, Total Newton Iters: %d, Residual: %6.4e\n',a1,a2,r.norm());

			fprintf('%10s  %10s  %10s  %10s  %10s\n','Iter','Step Size','|rz|','|rl|','|rv|');
		end
	end

	function print_detailed_line(o,iter,step_len,r)
		if o.display_level == 3
			[rz,rl,rv] = r.norms();
			a1 = int32(iter);
			fprintf('%10d  %10e  %10e  %10e  %10e\n',a1,step_len,rz,rl,rv);
		end

	end

	function print_detailed_footer(o,tol,r)
		if o.display_level == 3
			fprintf('Exiting inner loop. Inner residual: %6.4e, Inner tolerance: %6.4e\n',r.norm(),tol);
		end
	end

	% update options
	function update_options(o,opts)
		if isfield(opts,'sigma')
			o.sigma = opts.sigma;
		end

		if isfield(opts,'itol_red_factor')
			o.itol_red_factor = opts.itol_red_factor;
		end

		if isfield(opts,'max_newton_iters')
			o.max_newton_iters = opts.max_newton_iters;
		end

		if isfield(opts,'max_prox_iters')
			o.max_prox_iters = opts.max_prox_iters;
		end

		if isfield(opts,'tol')
			o.tol = opts.tol;
		end

		if isfield(opts,'rtol')
			o.rtol = opts.rtol;
		end

		if isfield(opts,'dtol')
			o.dtol = opts.dtol;
		end

		if isfield(opts,'inf_tol')
			o.inf_tol = opts.inf_tol;
		end

		if isfield(opts,'check_infeasibility')
			o.check_infeasibility = opts.check_infeasibility;
		end

		if isfield(opts,'alpha')
			o.alpha = opts.alpha;
		end

		if isfield(opts,'beta')
			o.beta = opts.beta;
		end

		if isfield(opts,'eta')
			o.eta = opts.eta;
		end

		if isfield(opts,'lsmax')
			o.lsmax = opts.lsmax;
		end

		if isfield(opts,'max_inner_iters')
			o.max_inner_iters = opts.max_inner_iters;
		end

		if isfield(opts,'itol_max')
			o.itol_max = opts.itol_max;
		end

		if isfield(opts,'itol_min')
			o.itol_min = opts.itol_min;
		end

		if isfield(opts,'display_level')
			o.display_level = opts.display_level;
		end
	end
end


end