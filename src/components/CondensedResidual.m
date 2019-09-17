% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% Computes and stores residuals for QPs
% without equality constraints.
classdef CondensedResidual < handle
properties(Access = public)
  data;
  alpha;
  rz; % stationarity residual
  rv; % complimentarity residual

  % cached norms
  znorm = 0;
  lnorm = 0;
  vnorm = 0;
end

methods(Access = public)
	function o = CondensedResidual(data)
		o.data = data;
		[nz,~,nv] = data.ProblemSize();
		o.rz = zeros(nz,1);
		o.rv = zeros(nv,1);
		o.alpha = 0.95;
	end

	% Computes the residuals for the inner PFB problem
	function InnerResidual(o,x,xbar,sigma)
		o.rz = o.data.H(x.z) + o.data.f + o.data.AT(x.v) +...
		 sigma*(x.z-xbar.z);
		ys = x.y + sigma*(x.v - xbar.v);
		o.rv = o.PFB(ys,x.v);

	    o.znorm = norm(o.rz);
	    o.vnorm = norm(o.rv);
	end

	% Computes the norm of the natural residual
	function NaturalResidual(o,x)
		o.rz = ...
		o.data.H(x.z) + o.data.f + o.data.AT(x.v);
		o.rv = o.alpha*min(x.y,x.v) + (1-o.alpha)*max(0,x.y).*max(0,x.v);

	    o.znorm = norm(o.rz);
	    o.vnorm = norm(o.rv);
	end
	
	% Compute the merit function
	function m = Merit(o)
		m = 1/2*norm(o)^2;
	end

	% overload norm
	function x = norm(o)
		x = sqrt(o.znorm^2 + o.vnorm^2);
	end

	function [rz,rl,rv] = GetNorms(o)
		rz = o.znorm;
		rv = o.vnorm;
		rl = 0;
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
	end
end % public methods

end