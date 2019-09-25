% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

% This class computes Newton steps directions
% using direct sprase linear solvers.
classdef SparseDirectSolver < handle
  properties(Access = public)
    ldl_piv_tol = 0.001;
    ztol = 100*eps;
    alpha = 0.95;

    % Pick the solver mode
    solver_mode = 2;
    use_ordering = true;

  end % public properties

  properties(Access = private)
    % Workspaces for storing matrix factorizations
    L; % lower triangular factor
    U; % upper triangular factor
    D; % block diagonal portion
    S; % scaling vector
   
    p; % row permutation vector
    q; % column permutation vector

    data; % problem data
    p_amd; % permutation from a minimum degree ordering

    gamma;
    mus;

    sigma;
  end

  methods(Access = public)

    function o = SparseDirectSolver(data,opts)
      if nargin < 2
        opts = struct();
      end

      if isfield(opts,'solver_mode')
        o.solver_mode = opts.solver_mode;
      end

      if isfield(opts,'use_ordering')
        o.use_ordering = opts.use_ordering;
      end

      o.data = data;
      [nz,nl,nv] = ProblemSize(o.data);

      % Perform a structural analysis to compute
      % a minimum degree ordering before factorization
      % TODO(dliaomcp@umich.edu) Add for mode 1
      if(o.solver_mode == 3 && o.use_ordering)
        K = o.data.H_ + speye(nz) + ...
        o.data.A_'*o.data.A_ + o.data.G_'*o.data.G_;
        o.p_amd = amd(K);
      elseif(o.solver_mode == 2 && o.use_ordering)
        E = o.data.H_ + speye(nz) + o.data.A_'*o.data.A_;
        K = [E,o.data.G_';o.data.G_,-speye(nl)];
        o.p_amd = amd(K);
      else
        o.p_amd = 1:nv;
      end
    end

    function Factor(o,x,xbar,sigma)
      [nz,nl,nv] = ProblemSize(o.data);

      o.sigma = sigma;
      % Compute the barrier terms
      ys = x.y + sigma*(x.v - xbar.v);
      [o.gamma,mu] = o.dphi(ys,x.v);
      o.mus = sigma*o.gamma + mu;

      Gamma = o.gamma./o.mus;

      switch o.solver_mode
      case 1
        o.Mode1Factor(o.gamma,o.mus,sigma);
      case 2
        o.Mode2Factor(Gamma,sigma);
      case 3
        o.Mode3Factor(Gamma,sigma);
      otherwise
        o.Mode2Factor(Gamma,sigma);
      end
     
    end

    function Solve(o,r,dx)
      switch o.solver_mode
      case 1
        o.Mode1Solve(r,dx);
      case 2
        o.Mode2Solve(r,dx);
      case 3
        o.Mode3Solve(r,dx);
      otherwise
        o.Mode2Factor(r,dx);
      end
    end

  end % public methods

  methods(Access = private)
    function Mode1Factor(o,gamma,mus,sigma)
      % Factor the sparse matrix
      % [Hs  G' A']
      % [-G  S  0 ]
      % [-CA 0  D ]
      [nz,nl,nv] = ProblemSize(o.data);

      Hs = o.data.H_ + sigma*speye(nz);
      C = spdiags(gamma,0,nv,nv);
      D = spdiags(mus,0,nv,nv);
      K = [Hs,o.data.G_',o.data.A_';
      -o.data.G_,sigma*speye(nl),sparse(nl,nv);
      -C*o.data.A_,sparse(nv,nl),D];

      % TODO(dliaomcp@umich.edu) Add preordering here using AMD or something
      [o.L,o.U,o.p,o.q,o.S] = lu(K,'vector');
    end

    function Mode1Solve(o,r,dx)
      [nz,nl,nv] = ProblemSize(o.data);

      rhs = [r.rz;r.rl;r.rv];
      x = zeros(nz+nl+nv,1);

      rhs = o.S\rhs;
      x(o.q) = o.U\(o.L\rhs(o.p));

      dx.z = x(1:nz);
      dx.l = x(nz+1:nz+nl);
      dx.v = x(nz+nl+1:end);
      dx.y = o.data.b - o.data.A(dx.z);

    end

    function Mode2Factor(o,Gamma,sigma)
      [nz,nl,nv] = ProblemSize(o.data);
      % Factor the sparse matrix
      % [E  G']
      % [G -sI]
     
      Gamma = spdiags(Gamma,0,nv,nv); % make it a sparse diagonal matrix
      E = o.data.H_ + sigma*speye(nz) + o.data.A_'*Gamma*o.data.A_;
      K = [E,o.data.G_';o.data.G_,-sigma*speye(nl)];

      [o.U,o.D,o.p,o.S] = ...
      ldl(K(o.p_amd,o.p_amd),o.ldl_piv_tol,'upper','vector');
    end

    function Mode2Solve(o,r,dx)
      [nz,nl,nv] = ProblemSize(o.data);

      % Compute the reduced residuals.
      rhs = zeros(nz+nl,1);
      rhs(1:nz) = r.rz - o.data.AT(r.rv./o.mus);
      rhs(nz+1:end) = -r.rl;

      rhs = rhs(o.p_amd);
      % Solve the system.
      rhs = o.S*rhs;
      x = zeros(nz+nl,1);
      x(o.p) = o.U\(o.D\(o.U'\(rhs(o.p))));
      x = o.S*x;

      % Invert the AMD permutation
      pinv(o.p_amd) = 1:length(o.p_amd);
      x = x(pinv);

      dx.z = x(1:nz);
      dx.l = x(nz+1:end);

      % Recover inequality duals.
      dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus;
      % dy for linesearch
      dx.y = o.data.b - o.data.A(dx.z);
    end


    function Mode3Factor(o,Gamma,sigma)
      [nz,nl,nv] = ProblemSize(o.data);

      Gamma = spdiags(Gamma,0,nv,nv); % make it a sparse diagonal matrix
      B = o.data.H_ + sigma*speye(nz) + ...
      o.data.A_'*Gamma*o.data.A_ + o.data.G_'*o.data.G_/sigma;

      [o.U,~,o.p] = chol(B(o.p_amd,o.p_amd),'vector');
    end

    function Mode3Solve(o,r,dx)
      [nz,nl,nv] = ProblemSize(o.data);

      % Compute the reduced residuals.
      rhs = zeros(nz,1);
      rhs = r.rz - o.data.AT(r.rv./o.mus) - o.data.GT(r.rl/o.sigma);
      rhs = rhs(o.p_amd);
      % Solve.
      dx.z = zeros(nz,1);
      dx.z(o.p) = o.U\(o.U'\(rhs(o.p)));

      % Invert the AMD permutation
      pinv(o.p_amd) = 1:length(o.p_amd);
      dx.z = dx.z(pinv);

      % Recover dual step directions.
      dx.l = 1/o.sigma*(r.rl + o.data.G(dx.z));
      dx.v = (r.rv + o.gamma.*o.data.A(dx.z))./o.mus; 
      % dy for linesearch
      dx.y = o.data.b - o.data.A(dx.z);
    end

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