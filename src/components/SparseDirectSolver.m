% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% This class computes Newton steps directions
% using direct sprase linear solvers.
classdef SparseDirectSolver < handle
  properties(Access = public)

    % The outputs from MATLABs ldl function
    U; % upper triangular factor
    D; % block diagonal portion
    p; % permutation vector
    S; % scaling vector
    K;

    data; % problem data
    p_amd; % permutation from a minimum degree ordering

    gamma;
    mus;

    ldl_piv_tol = 0.001;
    ztol = 100*eps;
    alpha = 0.95;

  end % public properties

  methods(Access = public)

    function o = SparseDirectSolver(data)
      o.data = data;

      % TODO(dliaomcp@umich.edu)
      % Perform a structural analysis to compute
      % a minimum degree ordering before factorization?
      % https://www.mathworks.com/help/matlab/math/sparse-matrix-reordering.html
    end
    
    % function PerformStructuralAnalysis(o)
    %   % TODO(dliaomcp@umich.edu)
    %   % Perform a structural analysis to decide if
    %   % the inequality blocks should be eliminated?
    %   % Decide based on predicted fill in.
    % end

    function Factor(o,x,xbar,sigma)
      [nz,nl,nv] = ProblemSize(o.data);

      % Compute the barrier terms
      ys = x.y + sigma*(x.v - xbar.v);
      [o.gamma,mu] = o.dphi(ys,x.v);
      o.mus = sigma*o.gamma + mu;

      Gamma = o.gamma./o.mus;
      Gamma = spdiags(Gamma,0,nv,nv); % make it a sparse diagonal matrix

      % TODO(dliaomcp@umich.edu) Add switch for matrix structure here.

      % Build the sparse matrix
      % [E  G']
      % [G -sI]

      E = o.data.H_ + sigma*speye(nz) + o.data.A_'*Gamma*o.data.A_;
      K = [E,o.data.G_';o.data.G_,-sigma*speye(nl)];
      % TODO(dliaomcp@umich.edu) Apply matrix reordering here
      % if helpful.
      [o.U,o.D,o.p] = ldl(K,o.ldl_piv_tol,'upper','vector');
    end

    function Solve(o,r,dx)
      [nz,nl,nv] = ProblemSize(o.data);

      % Compute the reduced residuals.
      rhs = zeros(nz+nl,1);
      rhs(1:nz) = r.rz - o.data.AT(r.rv./o.mus);
      rhs(nz+1:end) = -r.rl;

      % Solve the system.
      x = zeros(nz+nl,1);
      x(o.p) = o.U\(o.D\(o.U'\(rhs(o.p))));

      dx.z = x(1:nz);
      dx.l = x(nz+1:end);

      % Recover inequality duals.
      dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus;
      % dy for linesearch
      dx.y = o.data.b - o.data.A(dx.z);
    end

  end % public methods

  methods(Access = private)
    % compute an element of the C-differential 
    function [gamma,mu] = dphi(o,a,b,nv)
      [~,~,q] = ProblemSize(o.data);
      ztol = o.ztol;
      % computes an element from the C differential
      r = sqrt(a.^2 + b.^2);
      S0 =  r <= ztol; % zero elements
      S1 = a > 0 & b> 0; % positive orthant
      S2 = ~(S1 | S0); % the rest

      gamma = zeros(q,1);
      mu = zeros(q,1);
      for i = 1:q
        if S0(i)
          d = 1/sqrt(2);
          gamma(i) = o.alpha*(1- d);
          mu(i) = o.alpha*(1- d);
        elseif S1(i)
          gamma(i) = o.alpha*(1- a(i)/r(i)) + (1-o.alpha)*b(i);
          mu(i) = o.alpha*(1- b(i)/r(i)) + (1-o.alpha)*a(i);
        else % S2
          gamma(i) = o.alpha*(1- a(i)/r(i));
          mu(i) = o.alpha*(1- b(i)/r(i));
        end
      end
    end
  end % private methods


end % class