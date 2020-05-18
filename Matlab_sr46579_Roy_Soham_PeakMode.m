clc;
clear;
close all;

% Given
Pout = 6.6e3;
V_ac = 230;
f = 60;
Vc = 400;
Iout = Pout/Vc;
delta_Vc_pp = 20;
delta_Vc  = delta_Vc_pp/2;

Vg = V_ac*sqrt(2);
D = 1-(Vg/Vc);
D_prime = 1-D;
Ig = Iout/D_prime;

P_L = 0.5*0.01*Pout;

perc_ripple = 0.2;

fs = 100e3; % We are still in safe operating area
Ts = 1/fs;

omega = 2*pi*f;
C = Pout/(omega*Vc*2*delta_Vc);

L = V_ac^2*(1-(sqrt(2)*V_ac/Vc))/(0.2*fs*Pout);
R = Vc^2/Pout;

%% Preliminary design

% Assume
V_M = 1;
Rs = 2;
H = 5/400;
k_x_new = 0.045;

s = tf('s');
M_a = (Vc-Vg)/L;
M1 = Vg/L;
M2 = M_a;
V_a = Rs*M_a*Ts;
F_m = fs/(M_a + (M1-M2)/2);
F_v = D*D_prime*Ts/(2*L);
G_c0 = F_m*Vc/((D'*Rs)*(1+(2*F_m*Vc/(D_prime^2*R))+(F_m*F_v*Vc/D_prime)));
omega_c = (D_prime/sqrt(L*C))*sqrt(1+(2*F_m*Vc)/(D_prime^2*R)+(F_m*F_v*Vc/D_prime));
Qc = D_prime*R*sqrt(C/L)*sqrt(1+(2*F_m*Vc)/(D_prime^2*R)+(F_m*F_v*Vc/D_prime))/(1+(R*C*F_m*Vc/L)-(F_m*F_v*Vc/D_prime));
G_vc_cpm = G_c0*(1-(s*L)/(D_prime^2*R))/(1+(s/(Qc*omega_c))+(s/omega_c)^2);
Tu_cpm = G_vc_cpm*k_x_new*H*H*Vg;

% From PID tuner
Kp_new = 26.8;
Ti_new = 0.00429;
G_cv_cpm = Kp_new*(1+(1/(s*Ti_new)))
T_cpm = Tu_cpm*G_cv_cpm;

%% Final design

Kp_latest = 42.3;
Ti_latest = 0.00535;
G_cv_cpm_new = Kp_latest*(1+(1/(s*Ti_latest)))
T_cpm_new = zpk(Tu_cpm*G_cv_cpm_new);

%% Plots

opts = bodeoptions;
opts.Grid = 'on';
opts.FreqUnits = 'Hz';
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 0;
opts.PhaseMatchingValue = 0;
opts.Xlim=[1 fs];

figure;
bodeplot(Tu_cpm,opts);
title('T_u (cpm)');

figure;
bodeplot(T_cpm,opts);
title('T (cpm)');

figure;
bodeplot(T_cpm_new,opts);
title('T_{new} (cpm)');