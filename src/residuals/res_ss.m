% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% a class to compute the residuals for the FBstab solver
% this one handles the single shooting case
classdef res_ss < handle


properties(Access = public)
	data;

	alpha;

	rz; % stationarity residual
	rv;
end


methods(Access = public)
	% constructor
	function o = res_ss(data)
		o.data = data;
		[nz,~,nv] = data.opt_sz();
		[nx,nu,nc,N] = data.sz();
		o.rz = zeros(nz,1);
		o.rv = zeros(nv,1);
		o.alpha = 0.95;
	end

	% computes the residuals for the inner PFB problem
	function calcres(o,x,xbar,sigma)
		o.rz = o.data.H(x.z) + o.data.f + o.data.AT(x.v) + sigma*(x.z-xbar.z);

		ys = x.y + sigma*(x.v - xbar.v);
		% penalized FB function
		o.rv = o.pfb(ys,x.v);
	end

	% computes the norm of the natural residual
	function nres(o,x)
		o.rz = ...
		o.data.H(x.z) + o.data.f + o.data.AT(x.v);
		o.rv = o.alpha*min(x.y,x.v) + (1-o.alpha)*max(0,x.y).*max(0,x.v);
	end
	
	% compute the merit function
	function m = merit(o)
		m = 1/2*norm(o)^2;
	end

	% overload norm
	function x = norm(o)
		x = norm(o.rz) + norm(o.rv);
	end

	function [rz,rl,rv] = norms(o)
		rz = norm(o.rz);
		rv = norm(o.rv);
		rl = 0;
	end

	function vcopy(o,x)
		o.rv = x.rv;
		o.rv = x.rv;
	end

	% PFB function
	function y = pfb(o,a,b)
		yfb = a + b - sqrt(a.^2 + b.^2);
		ypen = max(0,a).*max(0,b);

		y = o.alpha* yfb + (1-o.alpha)*ypen;
	end

	% overload negate
	function negate(o)
		o.rz = -o.rz;
		o.rv = -o.rv;
	end
end


end