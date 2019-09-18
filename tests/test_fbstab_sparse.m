clear all
close all
clc

mpc = load_servo_example(30);
qp = form_sparse_mpc_qp(mpc);

nz = length(qp.f);
nl = length(qp.h);
nv = length(qp.b);

x0.z = zeros(nz,1);
x0.l = zeros(nl,1);
x0.v = zeros(nv,1);

opts.display_level = 2;

[x,out] = fbstab_sparse(x0,qp,opts);
