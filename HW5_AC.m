% Kayla Wong
% AERO421 - HW5

close all; clear; clc

J=diag([1200 2000 2800]); %kg*m2
%% Question 1

%design FSF control law for 2% settling time of 100s, damping ratio of 0.65
%and thrustuers as input torque. get the gain matrices 

t_s2 = 100; %sec, settling time
zeta=0.65; %damping ratio

%solve for natural frequency
w_n=4.4/(zeta*t_s2);

%damped frequency
w_d=w_n*sqrt(1-zeta^2);

%solve for max overshoot
M_p=exp(-pi*zeta/(sqrt(1-zeta^2)))*100;

%peak time
t_p=pi/w_d;

%rise time
beta=atan(sqrt(1-zeta^2)/zeta);
t_r=(pi-beta)/w_d; 

%gains matricies
kp=diag(2*w_n^2*diag(J));
kd=diag(2*zeta*w_n*diag(J));

disp("Proportional Gain matrix: ")
disp(kp)
disp("Derivative matrix: ")
disp(kd)

%% Question 2

%model in Simulink. model linear and non linear. does the non linear meet
%the performance requirements?

tspan = 120;

%initial conditions
epsilon_0=[0.2;-0.5;0.3];
eta_0=sqrt(1-epsilon_0.'*epsilon_0);
quat_0=[epsilon_0;eta_0];
w_0=[0.1; -0.05; 0.05]; %rad/s


% Case 1
quat_c=[0;0;0;1];
sim1=sim("HW5.slx");
q1_l=sim1.q_b_ECI;
w1_l=sim1.omega;
q1_nl=sim1.q_b_ECI_NL;
w1_nl=sim1.omega_NL;

figure(1)
subplot(2,1,1)
plot(q1_l.time, squeeze(q1_l.signals.values),'--',q1_nl.time, squeeze(q1_nl.signals.values),'-')
xlabel('Time (seconds)')
ylabel('Quaternion Components')
legend('\epsilon_x Linear','\epsilon_y Linear', '\epsilon_z Linear','\eta Linear', ...
'\epsilon_x Nonlinear', '\epsilon_y Nonlinear','\epsilon_z Nonlinear', '\eta Nonlinear')
grid on

subplot(2,1,2)
plot(w1_l.time, squeeze(w1_l.signals.values),'--',w1_nl.time, squeeze(w1_nl.signals.values),'-')
xlabel('Time (seconds)')
ylabel('Body Rates (rad/sec)')
legend('\omega_x Linear', '\omega_y Linear', '\omega_z Linear',...
    '\omega_x Nonlinear', '\omega_y Nonlinear', '\omega_z Nonlinear')
grid on
sgtitle("Nonlinear (solid) vs. Linearized (dashed) Response - Case 1")


% Case 2
epsilon_c2=[-0.2;0.4;0.2];
quat_c=[epsilon_c2;sqrt(1-epsilon_c2.'*epsilon_c2)];
sim2=sim("HW5.slx");
q2_l=sim2.q_b_ECI;
w2_l=sim2.omega;
q2_nl=sim2.q_b_ECI_NL;
w2_nl=sim2.omega_NL;

figure(2)
subplot(2,1,1)
plot(q2_l.time, squeeze(q2_l.signals.values),'--',q2_nl.time, squeeze(q2_nl.signals.values),'-')
xlabel('Time (seconds)')
ylabel('Quaternion Components')
legend('\epsilon_x Linear','\epsilon_y Linear', '\epsilon_z Linear','\eta Linear', ...
'\epsilon_x Nonlinear', '\epsilon_y Nonlinear','\epsilon_z Nonlinear', '\eta Nonlinear')
grid on

subplot(2,1,2)
plot(w2_l.time, squeeze(w2_l.signals.values),'-',w2_nl.time, squeeze(w2_nl.signals.values),'--')
xlabel('Time (seconds)')
ylabel('Body Rates (rad/sec)')
legend('\omega_x Linear', '\omega_y Linear', '\omega_z Linear',...
    '\omega_x Nonlinear', '\omega_y Nonlinear', '\omega_z Nonlinear')
grid on
sgtitle("Nonlinear (solid) vs. Linearized (dashed) Response - Case 2")


%% Question 3

%plot the torque required for both cases. assume thrusters are 2m from the
%COM of the sc and apply a pure moment about the desired axis of rotation.
%find a commercially available thruster that can generate the required
%thrust

d=2; %m


%Case 1
mc1_l=sim1.M_c;
mc1_nl=sim1.M_c_NL;

figure(3)
subplot(2,1,1)
plot(mc1_l.time, squeeze(mc1_l.signals.values).','--',mc1_nl.time,squeeze(mc1_nl.signals.values).','-')
ylabel('Torque (N-m)')
legend('M_{c,x} Linear','M_{c,y} Linear', 'M_{c,z} Linear', ...
'M_{c,x} Nonlinear', 'M_{c,y} Nonlinear','M_{c,z} Nonlinear')
grid on
title('Required Torque Commands - Case 1')

subplot(2,1,2)
plot(mc1_l.time, vecnorm(squeeze(mc1_l.signals.values).',2,2),'--', mc1_nl.time,...
    vecnorm(squeeze(mc1_nl.signals.values).',2,2))
xlabel('Time (seconds)')
ylabel('Norm of T_c (N-m)')
legend('||T|| Linear','||T|| Nonlinear')
grid on


%Case 2
mc2_l=sim2.M_c;
mc2_nl=sim2.M_c_NL;

figure(4)
subplot(2,1,1)
plot(mc2_l.time, squeeze(mc2_l.signals.values).','--',mc2_nl.time, squeeze(mc2_nl.signals.values).','-')
ylabel('Torque (N-m)')
legend('M_{c,x} Linear','M_{c,y} Linear', 'M_{c,z} Linear', ...
'M_{c,x} Nonlinear', 'M_{c,y} Nonlinear','M_{c,z} Nonlinear')
grid on
title('Required Torque Commands - Case 2')

subplot(2,1,2)
plot(mc2_l.time, vecnorm(squeeze(mc2_l.signals.values).',2,2),'--', mc2_nl.time,...
    vecnorm(squeeze(mc2_nl.signals.values).',2,2))
xlabel('Time (seconds)')
ylabel('Norm of T_c (N-m)')
legend('||T|| Linear','||T|| Nonlinear')
grid on



%Find required thrust from torque

%moment: T = F*d
% thrust req: F = T/d

%get max moment and force
Mmax_1_l=max(abs(squeeze(mc1_l.signals.values)), [],1);
Fmax_1_l=max(Mmax_1_l/d);

Mmax_1_nl=max(abs(squeeze(mc1_nl.signals.values)), [],1);
Fmax_1_nl=max(Mmax_1_nl/d);

Mmax_2_l=max(abs(squeeze(mc2_l.signals.values)), [],1);
Fmax_2_l=max(Mmax_2_l/d);

Mmax_2_nl=max(abs(squeeze(mc2_nl.signals.values)), [],1);
Fmax_2_nl=max(Mmax_2_nl/d);

disp("Case 1 Max Thrust Req: " + Fmax_1_nl + " N")
disp("Case 2 Max Thrust Req: " + Fmax_1_l + " N")

%Thruster: minimum 13 N, for 100second impulse
disp("Two fitting thrusters that meet the maximum thrust force required of" + ...
    "12.65 N, and a max estimated impulse time to meet the 100s settling " + ...
    "time are the HT-25N-SP by RAFAEL (9.5-28 N, 205-220s impulse), or" + ...
    "the Turn-Key Propulsion System by Benchmark Space Systems (10N, " + ...
    "150-170s impulse.")

