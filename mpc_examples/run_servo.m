% This file is part of the fbstab-matlab library
%
% https://github.com/dliaomcp/fbstab-matlab
%
% and is subject to the BSD-3-Clause license 


close all
clear all
clc



[mpc,sys,simin] = load_servo_example(25);



[t,U,X,Y] = sim_lti_mpc(sys,simin,mpc);

ytrg = simin.ytrg;
ymax = simin.ymax;
umax = simin.umax;

figure();
o = ones(size(t));
ax(1) = subplot(3,1,1);
plot(t,rad2deg(Y(1,:)),t,o*rad2deg(ytrg),'--');
legend('signal','target');
ylabel('Shaft Position [\circ]');
ax(2) = subplot(3,1,2);
plot(t,Y(2,:),t,o*ymax,'k--',t,-o*ymax,'k--');
ylabel('Shaft Torque [Nm]');
ax(3) = subplot(3,1,3);
h = stairs(t,U);
hold on 
plot(t,umax*o,'k--',t,-umax*o,'k--');
for i =1:length(h)
	h(i).LineWidth = 1.5;
end
ylabel('Input Voltage [V]')
linkaxes(ax,'x');
xlim([0,simin.tfinal]);
xlabel('Time [s]');