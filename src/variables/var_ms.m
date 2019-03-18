% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% class to hold a primal-dual variable
classdef var_ms < handle

properties(Access = public)
	z; % primal
	l; % equlity duals
	v; % inequlity duals
	y; % inequality margin

	% data object
	data;

end

methods(Access = public)
	
	function o = var_ms(data)
		o.data = data;
		o.zeros();
	end

	% loads a pdvar from a struct
	function scopy(o,s)
		% initialize from a struct
		o.z = s.z;
		o.l = s.l;
		o.v = s.v;
		o.y = o.data.b - o.data.A(s.z);
	end

	% write to a struct
	function s = swrite(o)
		s = struct();
		s.z = o.z;
		s.l = o.l;
		s.v = o.v;
	end

	% return a pdvar with all zeros
	function zeros(o)
		[nz,nl,nv] = o.data.opt_sz();
		o.z = zeros(nz,1);
		o.l = zeros(nl,1);
		o.v = zeros(nv,1);
		o.y = o.data.b;
	end

	% return a pdvar filled with all ones
	function ones(o)
		[nz,nl,nv] = o.data.opt_sz();
		o.z = ones(nz,1);
		o.l = ones(nl,1);
		o.v = ones(nv,1);
		o.y = o.data.b - o.data.A(o.z);
	end

	% performs a value copy i.e.,
	% copies x <- y
	function vcopy(x,y)
		x.z = y.z;
		x.l = y.l;
		x.v = y.v;
		x.y = y.y;
	end

	% compute y <- y + ax
	% for a scalar
	function axpy(y,x,a)
		y.z = y.z + a*x.z;
		y.l = y.l + a*x.l;
		y.v = y.v + a*x.v;
		y.y = y.y + a*(x.y - y.data.b);
	end

	% projects the dual variable onto the nonnegative orthant
	function dual_proj(o)
		o.v = max(0,o.v);
	end

	% overload norm
	function x = norm(o,ns)
		if nargin > 1
			x = norm(o.z,ns) + norm(o.l,ns) + norm(o.v,ns);
		else
			x = norm(o.z) + norm(o.l) + norm(o.v);
		end
	end

end

end