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
%
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