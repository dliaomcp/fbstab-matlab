% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% a class to compute the residuals for the FBstab solver
% this one handles the multiple shooting case
classdef FullResidual < handle
properties(Access = public)
 rz; % stationarity residual
 rl; % equality residual
 rv; % complimentarity residual
 data;
 alpha;

 % cached norms
 znorm = 0;
 lnorm = 0;
 vnorm = 0;
end

methods(Access = public)
	function o = FullResidual(data)
		o.data = data;
		[nz,nl,nv] = data.ProblemSize();
		o.rz = zeros(nz,1);
		o.rl = zeros(nl,1);
		o.rv = zeros(nv,1);
		o.alpha = 0.95;
	end

	% Computes the residuals for the inner PFB problem.
	function InnerResidual(o,x,xbar,sigma)
		o.rz = o.data.H(x.z) + o.data.f ...
		+ o.data.GT(x.l) + o.data.AT(x.v) + sigma*(x.z - xbar.z);
		o.rl = o.data.h - o.data.G(x.z) + sigma*(x.l - xbar.l);
		ys = x.y + sigma*(x.v - xbar.v);
		o.rv = o.PFB(ys,x.v);

		o.znorm = norm(o.rz);
		o.lnorm = norm(o.rl);
		o.vnorm = norm(o.rv);
	end

	% Computes the norm of the natural residual.
	function NaturalResidual(o,x)
		o.rz = o.data.H(x.z) + o.data.f ...
		+ o.data.GT(x.l) + o.data.AT(x.v);
		o.rl = o.data.h - o.data.G(x.z);
		o.rv = o.alpha*min(x.y,x.v); + (1-o.alpha)*max(0,x.y).*max(0,x.v);

		o.znorm = norm(o.rz);
		o.lnorm = norm(o.rl);
		o.vnorm = norm(o.rv);
	end

	% Computes the merit function.
	function m = Merit(o)
		m = 1/2*norm(o)^2;
	end

	% Computes the norm of the residual vector.
	function y = norm(o)
		y = sqrt(o.znorm^2 + o.lnorm^2 + o.vnorm^2);
	end

	function [rz,rl,rv] = GetNorms(o)
		rz = o.znorm;
		rv = o.vnorm;
		rl = o.lnorm;
	end

	% PFB function
	function y = PFB(o,a,b)
		yfb = a + b - sqrt(a.^2 + b.^2);
		ypen = max(0,a).*max(0,b);

		y = o.alpha* yfb + (1-o.alpha)*ypen;
	end

	% overload negate
	function Negate(o)
		o.rz = -o.rz;
		o.rv = -o.rv;
		o.rl = -o.rl;
	end
end % methods

end