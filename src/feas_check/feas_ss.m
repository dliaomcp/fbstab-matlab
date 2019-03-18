% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


% performs infeasibility checking for multiple shooting QPs
classdef feas_ss < handle

 properties(Access = public)
 	z1;
 	v1;

 	data;
 end


 methods(Access = public)

 	function o = feas_ss(data)
 		o.data = data;

 		[nz,nl,nv] = data.opt_sz();
 		o.z1 = zeros(nz,1);
 		o.v1 = zeros(nv,1);

 	end

 	function [primal,dual] = check_feas(o,x,tol)
 		[nz,nl,nv] = o.data.opt_sz();
 		primal = true;
 		dual = true;

 		% dual feasibility
 		w = norm(x.z,'inf');

 		o.z1 = o.data.H(x.z);
 		o.v1 = o.data.A(x.z);

 		d1 = norm(o.z1);
 		d2 = dot(o.data.f,x.z);
 		d3 = max(o.v1);

 		if d1 <= tol*w && d2 < 0 && d3 <= 0 && w > eps
 			dual = false;
 		end

 		% primal feasibility
 		u = norm(x.v,'inf');

 		o.z1 = zeros(nz,1) + o.data.AT(x.v);
 		p1 = norm(o.z1,'inf');

 		p2 = dot(x.v,o.data.b);

 		if p1 <= tol*u && p2 < 0 && u > eps
 			primal = false;
 		end
 	end

 end

end