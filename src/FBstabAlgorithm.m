% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

classdef FBstabAlgorithm < handle
properties(Access = private)
	% algorithm parameters
	sigma = 1000*sqrt(eps);
	sigma_max = eps^(1/4);
	sigma_min = 1e-10;
	max_newton_iters = 500;
	max_prox_iters = 100;
	tol = 1e-6;
	rtol = 1e-12;
	dtol = 1e-12;
	inf_tol = 1e-8;
	check_infeasibility = true;
	alpha = 0.95;
	beta = 0.9;
	eta = 1e-8;
	lsmax = 40;
	max_inner_iters = 50;
	itol_max = 1e-1;
	itol_min = 1e-12;
	itol_red_factor = 1/5;
	use_nonmonotone_linesearch = false;
	itol_relaxed = sqrt(eps);

	% display settings
	% 0 = none
	% 1 = final
	% 2 = outer iters
	% 3 = inner iters
	display_level = 0;

	% constituent classes
	data; % qp data object
	linear_solver; % linear system object
	feas_checker; % feasibility checker
	% residual objects
	rk; 
	ri;
	% variable objects
	xi;
	xk;
	xp;
	dx;

end % properties

properties(Access = private)
	fnorm;
	hnorm;
	bnorm;

	ztol = 0;
	ltol = 0;
	vtol = 0;
end

methods(Access = public)
% constructor 
	function o = FBstabAlgorithm(data,linsol,feas,rk,ri,xk,xi,xp,dx,opts)
		if nargin < 5
			opts = struct();
		else
			o.UpdateOptions(opts);
		end

		o.data = data;
		o.linear_solver = linsol;
		o.feas_checker = feas;
		o.rk = rk;
		o.ri = ri;

		o.xk = xk;
		o.xi = xi;
		o.xp = xp;
		o.dx = dx;

		% harmonize alpha
		o.rk.alpha = o.alpha;
		o.ri.alpha = o.alpha;
		linear_solver.alpha = o.alpha;

		o.fnorm = norm(o.data.f,'inf');
		o.hnorm = norm(o.data.h,'inf');
		o.bnorm = norm(o.data.b,'inf');
	end

  % Solves the QP.
	function [x,out] = Solve(o,x0)
		[nz,nl,nv] = ProblemSize(o.data);

		%% initialization *************************************
		xk = o.xk; % Outer loop variable
		xi = o.xi; % Inner loop variable
		dx = o.dx; % Workspace
		xp = o.xp; % Workspace

		xk.StructRead(x0);
		dx.ones();

		rk = o.rk; % Outer loop
		ri = o.ri; % Inner loop
		rk.NaturalResidual(xk);

		% Get quality of the initial guess.
		E0 = norm(rk);
		dE = 1e6;

		% Initialize output structure.
		out = struct();
		out.eflag = -1;
		out.newton_iters = 0;
		out.prox_iters = 0;
		out.res = E0;

		% Initialize tolerance
		itol = o.saturate(E0,o.itol_min,o.itol_max);

		% initial regularization parameter
		sigma = o.sigma;
		
		o.PrintIterationHeader();

		% Begin the main proximal loop.
		for k = 1:o.max_prox_iters
			% Termination check.
			rk.NaturalResidual(xk);
			Ek = rk.norm();

			tol = o.GetTolerances(xk);
			[rz,rl,rv] = rk.GetNorms();
			if rz <= tol.z && rl <= tol.l && rv <= tol.v
			% if Ek <= o.tol + E0*o.rtol
				out.eflag = 0;
				break;
			elseif out.newton_iters > o.max_newton_iters
				out.eflag = -1;
				break;
			elseif dE <= o.dtol
				out.eflag = -2;
			end

			o.PrintDetailedHeader(out.prox_iters,out.newton_iters,rk,sigma);
			o.PrintIterationLine(out.prox_iters,out.newton_iters,rk,ri,itol,sigma);
			
			Ek_old = Ek;
			% Call inner solver *************************************
			[Ei,Ek,out,flag] = o.SolveSubproblemFBRS(sigma,itol,Ek,out);

			% Reduce sigma if the iteration is deemed to have stalled or is successfull.
			if(flag == 1 || flag == 2)
				sigma = max(sigma/10,o.sigma_min);
				itol = o.saturate(o.itol_red_factor*itol,o.itol_min,Ek);
			elseif(flag == 5)
				sigma = max(sigma/10,o.sigma_min);
			end
			% Increase sigma after an iteration timeout that indicates 
			% that the inner solver failed.
			if(flag == 3)
				sigma = o.saturate(sigma*10,o.sigma_min,o.sigma_max);
				itol = o.saturate(itol/o.itol_red_factor,o.itol_min,Ek);
			end
			if(flag == -1)
				break;
			end

			% Accept the update if it decreases the problem residual
			if Ek <= Ek_old
				% dx = xi - xk	
				dx.Copy(xi);
				dx.axpy(-1,xk);
				xk.Copy(xi);
			end
			% Feasibility Check
			if o.check_infeasibility && Ek <= Ek_old
				[primal,dual] = o.feas_checker.CheckFeasibility(dx,o.inf_tol);
				if ~primal && ~dual
					out.eflag = -5;
					x = dx.StructWrite();
					out.res = rk.norm();
					o.PrintFinal(out,rk);
					return;
				elseif ~primal
					out.eflag = -3;
					x = dx.StructWrite();
					out.res = rk.norm();
					o.PrintFinal(out,rk);
					return;
				elseif ~dual
					out.eflag = -4;
					x = dx.StructWrite();
					out.res = rk.norm();
					o.PrintFinal(out,rk);
					return;
				end
			end
			out.prox_iters = out.prox_iters +1;
		end % prox loop

		% output solution
		rk.NaturalResidual(xk);
		out.res = rk.norm();
		x = xk.StructWrite();
		o.PrintFinal(out,rk);
	end % end solve

	function [Ei,Ek,out,flag] = SolveSubproblemFBRS(o,sigma,itol,Eouter,out)
		x = o.xi;
		xbar = o.xk;
		dx = o.dx;
		xp = o.xp;
		ri = o.ri;
		rk = o.rk;

		% Different modes:
		% 1: Everything OK -> reduce tolerances
		% 2: Stall, subproblem solved but no progress is made -> reduce sigma
		% 3: Out of iterations: don't reduce itol
		% 4: itol_min hit -> usually happens during infeasibility detection
		% 5: having trouble hitting full target but can hit relaxed -> reduce sigma
		% -1: Out of newton iterations
		flag = 1;
		merit_buffer = zeros(5,1);
		t = 1;
		Ei = 0;

		x.Copy(xbar);
		for j = 1:o.max_inner_iters
			% Convergence check.
			ri.InnerResidual(x,xbar,sigma);
			rk.NaturalResidual(x);

			Ei = ri.norm();
			Ek = rk.norm();
			dx.Copy(x);
			dx.axpy(-1,xbar);

			if (Ei <= itol*min(1,norm(dx)) && norm(rk) >= Eouter)
				flag = 2;
				o.PrintDetailedFooter(itol,ri);
				return;
			elseif (Ei <= itol*min(1,norm(dx)) && norm(rk) < Eouter)
				flag = 1;
				o.PrintDetailedFooter(itol,ri);
				return;
			elseif (Ei <= o.itol_min*100 && norm(rk) <= Eouter)
				o.PrintDetailedFooter(itol,ri);
				flag = 5;
				return;
			elseif (Ei <= o.itol_min)
				o.PrintDetailedFooter(itol,ri);
				flag = 4;
				return;
			end
			o.PrintDetailedLine(j,t,ri);

			if out.newton_iters >= o.max_newton_iters
				out.eflag = -1;
				rk.NaturalResidual(xbar);
				out.res = rk.norm();
				x = xbar.StructWrite();
				o.PrintFinal(out,rk);
				flag = -1;
				return
			end

			% form and factor the iteration matrix
			o.linear_solver.Factor(x,xbar,sigma);

			% Solve the system for the rhs -ri and  store the result in dx.
			ri.Negate();
			o.linear_solver.Solve(ri,dx);

			% Linesearch
			merit_buffer = circshift(merit_buffer,1);
			merit_buffer(1) = 1/2*Ei^2;
			if(o.use_nonmonotone_linesearch)
				m0 = max(merit_buffer);
			else
				m0 = merit_buffer(1);
			end

			t = 1;
			for i = 1:o.lsmax
				% compute xp = x + t*dx
				xp.Copy(x);
				xp.axpy(t,dx);
				ri.InnerResidual(xp,xbar,sigma);
				mp = 1/2*norm(ri)^2;
				if mp <= m0 - 2*t*o.eta*merit_buffer(1)
					break;
				else
					t = o.beta*t;
				end
			end
			x.axpy(t,dx);  % x = x + t*dx
			out.newton_iters = out.newton_iters+1;
		end % pfb loop
		% Iteration timeout flag
		flag = 3;
		x.ProjectDuals();
	end

	% Computes the tolerances at a given point.
	function tol = GetTolerances(o,x)
		o.ztol = o.fnorm + 1;
		o.ltol = o.hnorm + 1;
		o.vtol = o.bnorm + 1;

		o.ztol = o.rtol*o.ztol + o.tol;
		o.ltol = o.rtol*o.ltol + o.tol;
		o.vtol = o.rtol*o.vtol + o.tol;

		tol.z = o.ztol;
		tol.l = o.ltol;
		tol.v = o.vtol;
	end
	% printing
	function PrintFinal(o,out,r)
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

			[rz,rl,rv] = r.GetNorms();

			a1 = int32(out.prox_iters);
			a2 = int32(o.max_prox_iters);

			fprintf('Proximal iterations: %d out of %d\n', a1, a2);

			a3 = int32(out.newton_iters);
			a4 = int32(o.max_newton_iters);
			fprintf('Newton iterations: %d out of %d\n', a3, a4);
			fprintf('%10s  %10s  %10s  %10s  %10s  %10s\n ','|rz|','|rl|','|rv|','z tol','l tol','v tol');
			fprintf('%10.4e  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e\n',rz,rl,rv,o.ztol,o.ltol,o.vtol);
		end
	end

	function PrintIterationHeader(o)
		if o.display_level == 2
			fprintf('%12s %12s %12s %12s %12s %12s %12s %12s\n','prox iter','newton iters','|rz|','|rl|','|rv|','Inner res','Inner tol','sigma');
		end
	end

	function PrintIterationLine(o,prox_iters,newton_iters,rk,ri,itol,sigma)
		if o.display_level == 2
			[rz,rl,rv] = rk.GetNorms();
			a1 = int32(prox_iters);
			a2 = int32(newton_iters);
			fprintf('%12d %12d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n',a1,a2,rz,rl,rv,ri.norm(),itol,sigma);
		end
	end

	function PrintDetailedHeader(o,prox_iters, newton_iters, r,sigma)
		if o.display_level == 3
			a1 = int32(prox_iters);
			a2 = int32(newton_iters);
			fprintf('Begin Prox Iter: %d, Total Newton Iters: %d, Residual: %6.4e, sigma: %6.4e\n',a1,a2,r.norm(),sigma);

			fprintf('%10s  %10s  %10s  %10s  %10s\n','Iter','Step Size','|rz|','|rl|','|rv|');
		end
	end

	function PrintDetailedLine(o,iter,step_len,r)
		if o.display_level == 3
			[rz,rl,rv] = r.GetNorms();
			a1 = int32(iter);
			fprintf('%10d  %10e  %10e  %10e  %10e\n',a1,step_len,rz,rl,rv);
		end

	end

	function PrintDetailedFooter(o,tol,r)
		if o.display_level == 3
			fprintf('Exiting inner loop. Inner residual: %6.4e, Inner tolerance: %6.4e\n',r.norm(),tol);
		end
	end

	% update options
	function UpdateOptions(o,opts)
		if isfield(opts,'sigma')
			o.sigma = opts.sigma;
		end

		if isfield(opts,'sigma_min')
			o.sigma_min = opts.sigma_min;
		end

		if isfield(opts,'sigma_max')
			o.sigma_max = opts.sigma_max;
		end

		if isfield(opts,'use_nonmonotone_linesearch')
			o.use_nonmonotone_linesearch = opts.use_nonmonotone_linesearch;
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
end % methods

methods(Access = private)
	% Projects x onto [a,b]
	function y = saturate(o,x,a,b)
		y = min(x,b);
		y = max(y,a);
	end
end

end % classdef