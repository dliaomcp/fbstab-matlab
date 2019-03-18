% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


clear all
close all
clc

[mpc,sys,simin] = load_hcw_example(40);

[t,U,X,Y] = sim_lti_mpc(sys,simin,mpc);


umax = simin.umax;
vmax = simin.vmax;

figure();
subplot(3,1,1);
plot(t,X(1:3,:));
ylabel('Position [km]');
legend('Radial','Along Track','Across Track','orientation','horizontal');

o = ones(size(t));
subplot(3,1,2);
plot(t,X(4:6,:),t,vmax*o,'k--',t,-vmax*o,'k--');
ylabel('Velocity [km/s]');

subplot(3,1,3);
h = stairs(t,U');
for i =1:length(h)
	h(i).LineWidth = 1.5;
end
hold on
plot(t,umax*o,'k--',t,-umax*o,'k--');
ylabel('Impulse [km/s]');
xlabel('Time [s]');