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
k_x = Ig*Rs/Vg;

s = tf('s');

G_id0 = (2*Vc)/(D_prime^2*R);
f0 = D_prime/(2*pi*sqrt(L*C));
omega0 = 2*pi*f0;
f_zi = 1/(pi*R*C);
omega_zi = 2*pi*f_zi;
Q = D_prime*R*sqrt(C/L);

G_id = zpk(G_id0*(1+(s/omega_zi))/(1+(s/(Q*omega0))+(s/omega0)^2));

T_iu = zpk((1/V_M)*G_id*Rs);

f_ci = (1/6)*fs;
omega_ci = 2*pi*f_ci;
boost = 52; % desired phase margin
k = tan(deg2rad((boost/2) + 45));
omega_zci = omega_ci/k;
omega_pci = omega_ci*k;

G_cm = (L*omega_ci*V_M)/(Vc*Rs);
G_ci = G_cm*(1+(omega_zci/s))/(1+(s/omega_pci))

T_i = zpk(T_iu*G_ci);

opts = bodeoptions;
opts.Grid = 'on';
opts.FreqUnits = 'Hz';
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 0;
opts.PhaseMatchingValue = 0;
opts.Xlim=[1 fs];

% Voltage loop

G_vc0 = (D_prime*R)/(2*Rs);
omega_zv = (D_prime)^2*R/L;
f_zv = omega_zv/(2*pi);
omega_pv = 2/(R*C);
f_pv = omega_pv/(2*pi);

G_vc = G_vc0*(1-(s/omega_zv))/(1+(s/omega_pv));

T_vu = G_vc*H;

f_cv = 30;
omega_cv = 2*pi*f_cv;
omega_zcv = omega_cv/3;

G_vm = (2*pi*f_cv*C*Rs)/(D_prime*H);
G_cv = G_vm*(1+(omega_zcv/s))

T_v = zpk(T_vu*G_cv);

%% Final design
f_h = 120;
f_cv_new = 20;
omega_cv_new = 2*pi*f_cv_new;

% From PID tuner
Kp = 26.8;
Ti = 0.00429;
G_cv_new = Kp*(1+(1/(s*Ti)))
T_v_new = zpk(T_vu*G_cv_new);

%% Plots

figure;
bodeplot(T_iu,opts);
title('T_{iu}');

figure;
bodeplot(T_i,opts);
title('T_i');

figure;
bodeplot(T_vu,opts);
title('T_{vu}');

figure;
bodeplot(T_v,opts);
title('T_v');

figure;
bodeplot(T_v_new,opts);
title('T_{v,new}');

%% Asymptotes

figure;
asymp(T_iu, opts);
title('T_{iu}');