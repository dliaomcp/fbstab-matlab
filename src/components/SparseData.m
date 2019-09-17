% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% This class provides storage and implements methods for working with
% general sparse QPs.
classdef SparseData < handle

methods(Access = public)
 function o = SparseData(H,G,A,f,h,b)
 	%% obj = SparseData(H,G,A,f,h,b)
 	% Copies in problem data assuming
 	% H, G, and A are sparse matrices and
 	% f,h, and b are vectors.
 	o.nz = length(f);
 	o.nl = length(h);
 	o.nv = length(b);

 	o.H = H;
 	o.A = A;
 	o.G = G;
 	o.f = f;
 	o.h = h;
 	o.b = b;

 	% TODO(dliaomcp@umich.edu) add assertions that the sizes all match up.
 end

 function [nz,nl,nv] = ProblemSize(obj)
 	nz = obj.nz;
 	nl = obj.nl;
 	nv = obj.nv;
 end

 function y = H(obj,v)
 	y = obj.H*v;
 end

 function y = A(obj.v)
 	y = obj.A*v;
 end

 function y = G(obj.v)
 	y = obj.G*v;
 end

 function y = AT(obj.v)
 	y = (obj.A')*v;
 end

 function y = GT(obj.v)
 	y = (obj.G')*v;
 end
end

properties(Access = private)
 nz = 0;
 nl = 0;
 nv = 0;

 H;
 A;
 G;
 f;
 b;
 h;

end

end % class