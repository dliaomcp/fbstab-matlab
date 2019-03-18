close all
clear all
clc



[mpc,sys,simin] = load_copoly_example(70);



[t,U,X,Y] = sim_lti_mpc(sys,simin,mpc);

figure();
subplot(2,1,1)
plot(t,Y);
legend('y_1','y_2','y_3','y_4','orientation','horizontal');
ylabel('Normalized Outputs');

o = ones(size(t));
subplot(2,1,2);
plot(t,U);
legend('u_1','u_2','u_3','u_4','u_5','orientation','horizontal');
ylabel('Normalized Inputs');
xlabel('Time');