% @Author: Dominic Liao-McPherson
% @Date:   2018-10-30 13:20:48
% @Last Modified by:   Dominic Liao-McPherson
% @Last Modified time: 2019-04-25 17:15:40


% a function to call the multiple shooting version of fbstab
% using the banded LDL' factorization approach for the linear systems
function [x,u,out] = fbstab_ldl_sp(x0,mpc,opts)
	if nargin < 3
		opts = struct();
	end

	data = data_ms(mpc.Q,mpc.R,mpc.S,mpc.q,mpc.r,mpc.A,...
		mpc.B,mpc.c,mpc.x0,mpc.E,mpc.L,mpc.d);

	linsys = ldlt_sp(data);
	feas = feas_ms(data);
	r1 = res_ms(data);
	r2 = res_ms(data);

	x1 = var_ms(data);
	x2 = var_ms(data);
	x3 = var_ms(data);
	x4 = var_ms(data);
	% solver object
	fbstab = fbstab_algo(data,linsys,feas,r1,r2,x1,x2,x3,x4,opts);

	% call the solver
	[x,out] = fbstab.solve(x0);
	[nx,nu] = size(mpc.B(:,:,1));
	u = x.z(nx+1:nx+nu);

end