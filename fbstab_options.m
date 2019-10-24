%% fbstab_options: Generate default options for fbstab
% opts = fbstab_options(option_type) returns a structure of options where
% option_type is one of 'default', or 'reliable'.
% The default options are recommended when solving sequences of problems quickly,
% e.g., when solving model predictive control problems.
% The reliable options are designed to make the algorithm 
% more robust and tends to reduce the iteration counts for harder problems where
% accelerating heuristics tend to hurt rather than help.
% The list of options is as follows:
%
% display_level: controls printed output
%  0: Silent
%  1: Final 
%  2: Outer Iterations
%  3: All Iterations 
%
% sigma: Initial stabilization parameter
% sigma_min: Minimum stabilization parameter
% sigma_max: Maximum stabilization parameter
% max_newton_iters: Maximum number of Newton iterations allowed
% max_prox_iters: Maximum number of proximal point iterations
% tol: absolute tolerance
% rtol: relative tolerance
% dtol: stall tolerance
% inf_tol: infeasibility tolerance
% check_infeasibility: enable or disable infeasibility checking
% alpha: penalized Fischer-Burmeister parameter
% beta: backtracking linesearch parameter, i.e., step <- beta*step
% eta: sufficient decrease parameter
% lsmax: maximum number of linesearch iterations
% max_inter_iters: maximum number of newton iterations per proximal iterations
% itol_max: maximum subproblem tolerance
% itol_min: minimum subproblem tolerance
% itol_red_factor: reduce the itol by this factor if subproblem is successful
% use_nonmonotone_linesearch: enable or disable nonmonotone linesearch
% 
% fbstab_sparse only:
%   solver_mode: 
%     1: solve the asymmetric newton step system with sparse LU
%     2: solve the quasidefinite reduced newton system with LDL'
%     3: solve the normal equations with LL'
%   use_ordering: enable or disable a prefactorization amd ordering
%
% fbstab_mpc only:
%   linear_solver: Which linear solver to use
%     ric: Ricatti recursion based method
%     pcg: Implicit conjugate gradient
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

function opts = fbstab_options(option_type)
  if nargin < 1
    option_type = 'default';
  end
  opts.sigma = sqrt(eps);
  opts.sigma_max = opts.sigma;
  opts.sigma_min = 1e-10;
  opts.max_newton_iters = 200;
  opts.max_prox_iters = 30;
  opts.tol = 1e-6;
  opts.rtol = 1e-12;
  opts.dtol = 1e-12;
  opts.inf_tol = 1e-8;
  opts.check_infeasibility = true;
  opts.alpha = 0.95;
  opts.beta = 0.75;
  opts.eta = 1e-8;
  opts.lsmax = 30;
  opts.max_inner_iters = 50;
  opts.itol_max = 1e-1;
  opts.itol_min = 1e-12;
  opts.itol_red_factor = 1/5;
  opts.use_nonmonotone_linesearch = true;
  opts.use_ordering = true;
  opts.solver_mode = 3;
  opts.linear_solver = 'ric';
  opts.display_level = 2;

  if strcmp(option_type,'reliable')
    opts.beta = 0.9;
    opts.lsmax = 40;
    opts.sigma = 1000*sqrt(eps);
    opts.sigma_max = eps^(1/4);
    opts.sigma_min = 1e-10;
    opts.max_newton_iters = 500;
    opts.max_prox_iters = 100;
    opts.use_nonmonotone_linesearch = false;
    opts.solver_mode = 2;
  end
end