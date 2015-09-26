clc;
clear all;

m = 1.650; % mass = 0.650[kg]
g = 9.81; %[m/s^2]
Ix = 7.5e-3; % inertia on x axis [kg*m^2]
Iy = 7.5e-3; % inertia on y axis [kg*m^2]
Iz = 1.3e-2; % inertia on z axis [kg*m^2]
Inerita_Mat = diag([Ix, Iy, Iz]);

b = 3.13e-5; % thrust coefficient [N*s^2] 
d = 7.5e-7;  % drag coefficient [N*m*s^2]
Jr = 6e-5;   % rotor inertia    [kg*m^2]
L_arm = 0.23; % arm length      [m]

G_motor = 1;
%G_motor = tf([1],[0.01 1]);
%G_motor = tf([0.936],[0.178 1]);
Delay_Act = 0.0;


Rotor_State_Normal = diag([1 1 1 1]);                
Rotor_State_Fail = diag([1 1 1 0.5]);

Speed_Square_to_Thrust_N_Moment_Mat = [b b b b; 0 -b 0 b; -b 0 b 0; -d d -d d];

Sum_Rotor_Angular_Velocity_Mat = [-1 1 -1 1];

%%% Trajectory Tracking BackStepping Controller %%%
alpha1 = 10;
alpha2 = .3;
alpha3 = 1;
alpha4 = .1;
alpha5 = 1;
alpha6 = .1;

%%% Rotational Control PID L1 StateFB %%%

Am = - diag([10, 10, 10])*2;
Bm = diag([L_arm/Ix, L_arm/Iy, 1/Iz]);

Q = 1*eye(3);
P = lyap(Am',Q);

Ts = 0.005;
