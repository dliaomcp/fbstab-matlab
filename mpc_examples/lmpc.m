

function [u,s,res] = lmpc(x,xbar,ubar,qp,s)

	% set up the QP
	p.H = (qp.H + qp.H')/2;
	p.f = qp.h + qp.Tx*xbar + qp.Tu*ubar + qp.Tf*x;
	p.A = qp.G;
	p.b = qp.g + qp.Tb*x;


	[z,v,out] = fbstab(p,s.z,s.v);

	res = out.res;
	s.z = z;
	s.v = v;
	u = z(1:qp.nu);

end