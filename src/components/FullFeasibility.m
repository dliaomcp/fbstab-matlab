% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% Performs infeasibility checking for QPs with equality constraints.
classdef FullFeasibility < handle

properties(Access = private)
 	z1;
 	l1;
 	v1;

 	data;
end

methods(Access = public)
	function o = FullFeasibility(data)
		o.data = data;

		[nz,nl,nv] = data.ProblemSize();
		o.z1 = zeros(nz,1);
		o.l1 = zeros(nl,1);
		o.v1 = zeros(nv,1);

	end

	function [primal,dual] = CheckFeasibility(o,x,tol)

		primal = true;
		dual = true;

		% dual feasibility
		w = norm(x.z,'inf');

		o.z1 = o.data.H(x.z);
		o.l1 = o.data.G(x.z);
		o.v1 = o.data.A(x.z);

		d1 = norm(o.z1);
		d2 = dot(o.data.f,x.z);
		d3 = max(o.v1);
		d4 = norm(o.l1);

		if d1 <= tol*w && d2 < 0 && d3 <= 0 && d4 <= tol*w && w > eps
			dual = false;
		end

		% primal feasibility
		u = norm(x.l,'inf') + norm(x.v,'inf');

		o.z1 = o.data.GT(x.l) + o.data.AT(x.v);
		p1 = norm(o.z1,'inf');

		p2 = dot(x.v,o.data.b) + dot(x.l,o.data.h);

		if p1 <= tol*u && p2 < 0 && u > eps
			primal = false;
		end
	end

end

end