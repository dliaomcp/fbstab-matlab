% fbstab_sparse: Proximally Stabilized Semismooth Method for solving QPs
% [x,out] = fbstab_mpc(x0,qp,opts) attempts to solve the 
% following quadratic programming problem:
%
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az ≤ b
%
% Aside from convexity there are no assumptions made about the problem
% This method can detect unboundedness/infeasibility
% and can exploit arbitrary initial guesses. This implementation uses 
% dense linear algebra and dense factorizations.
% 
% Inputs: 
% x0:  A structure with the following fields
%   z: primal initial guess
%   l: equality dual initial guess
%   v: inequality dual initial guess
%
% qp: A structure with the following fields
%   H: n x n Sparse matrix
%   f: n x 1 
%   G: m x n Sparse matrix
%   h: m x 1
%   A: q x n Sparse matrix
%   b: q x 1
%
% opts: Options structure, see fbstab_options.m
% 
% Outputs:
% x: Structure with the following fields
%   z: primal solution
%   l: equality duals
%   v: inequality duals
%
% out: Structure with the following fields
%   prox_iters: number of proximal iterations 
%   newton_iters: number of newton iterations
%   res:  KKT residual
%   eflag: Exit flag
%      0: success
%     -1: Maximum number of iterations exceeded
%     -2: Algorithm stalled
%     -3: Problem is infeasible
%     -4: Problem is unbounded below (dual infeasible)
%     -5: Problem is primal and dual infeasible
%
% Details about the the FBstab algorithm can be found at:
% https://arxiv.org/abs/1901.04046
%
% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 

%Copyright 2018-2019 University of Michigan

% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [x,out] = fbstab_sparse(x0,qp,opts)
  if nargin < 3
      opts = struct();
  end
  % Remove infinities from the inequality constraints.
  n = length(qp.b);
  idx = ~isinf(qp.b);
  x0.v = x0.v(idx);

  data = SparseData(qp.H,qp.G,qp.A(idx,:),qp.f,qp.h,qp.b(idx));
  linsys = SparseDirectSolver(data,opts);
  feas = FullFeasibility(data);
  r1 = FullResidual(data);
  r2 = FullResidual(data);
  x1 = FullVariable(data);
  x2 = FullVariable(data);
  x3 = FullVariable(data);
  x4 = FullVariable(data);

  fbstab = FBstabAlgorithm(data,linsys,feas,r1,r2,x1,x2,x3,x4,opts);

  [x,out] = fbstab.Solve(x0);

  % Modify v by filling the indices corresponding to infinite inequalities 
  % with 0 (inactive).
  vnew = zeros(n,1);
  vnew(idx) = x.v;
  x.v = vnew;
end
