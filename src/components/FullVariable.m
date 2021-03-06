% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% This class implements primal-dual variables for QPs with 
% equality and inequality constraints.
classdef FullVariable < handle
properties(Access = public)
	z; % primal
	l; % equlity duals
	v; % inequlity duals
	y; % inequality margin
	data; % data object
end

methods(Access = public)
	function o = FullVariable(data)
		o.data = data;
		o.zeros();
	end

	% Loads a variable from a struct
	function StructRead(o,s)
		% initialize from a struct
		o.z = s.z;
		o.l = s.l;
		o.v = s.v;
		o.y = o.data.b - o.data.A(s.z);
	end

	% Write to a struct
	function s = StructWrite(o)
		s = struct();
		s.z = o.z;
		s.l = o.l;
		s.v = o.v;
	end

	% Fill with all zeros.
	function zeros(o)
		[nz,nl,nv] = o.data.ProblemSize();
		o.z = zeros(nz,1);
		o.l = zeros(nl,1);
		o.v = zeros(nv,1);
		o.y = o.data.b;
	end

	% Fill with all ones
	function ones(o)
		[nz,nl,nv] = o.data.ProblemSize();
		o.z = ones(nz,1);
		o.l = ones(nl,1);
		o.v = ones(nv,1);
		o.y = o.data.b - o.data.A(o.z);
	end

	% Performs a value copy i.e., x <- y
	function Copy(x,y)
		x.z = y.z;
		x.l = y.l;
		x.v = y.v;
		x.y = y.y;
	end

	% compute y <- y + ax
	% for a scalar
	function axpy(y,a,x)
		y.z = y.z + a*x.z;
		y.l = y.l + a*x.l;
		y.v = y.v + a*x.v;
		y.y = y.y + a*(x.y - y.data.b);
	end

	% projects the dual variable onto the nonnegative orthant
	function ProjectDuals(o)
		o.v = max(0,o.v);
	end

	% overload norm
	function x = norm(o)
		x = sqrt(norm(o.z)^2 + norm(o.l)^2 + norm(o.v)^2);
	end

end

end