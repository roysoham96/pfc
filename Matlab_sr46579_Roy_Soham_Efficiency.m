clear;
close all;
clc;

% Parameters
Ron = 40e-3; 
Coss = 85e-12;

% Full load to half load graph
Pout = [6.6e3 5.5e3 4.4e3 3.3e3];
V_ac = 230;
f = 60;
Vc = 400;
Iout = Pout./Vc;
delta_Vc_pp = 20;
delta_Vc  = delta_Vc_pp/2;

Vg = V_ac*sqrt(2);
D = 1-(Vg/Vc);
D_prime = 1-D;
Ig = Iout./D_prime;

P_L = Pout.*0.5*0.01;

fs = 100e3; % We are still in safe operating area
Ts = 1/fs;

eta_min = 0.978;

I_ac = Pout/(eta_min*V_ac);
I_switch_rms = I_ac*sqrt(1-(8*Vg)/(3*pi*Vc));
P_switch_con = I_switch_rms.^2.*Ron;

P_sw = 0.5*Coss*Vc^2*fs;

P_switch = P_switch_con + P_sw;

I_D_rms = I_ac*sqrt(8*Vg/(3*pi*Vc));
V_F_D1 = 1.5;
P_D = I_D_rms.*V_F_D1;

V_F_D2 = 0.795;
Iavg = 2*sqrt(2)*Pout./(pi*V_ac);
P_DBR = 2*Iavg.*V_F_D2;

P_losses = P_DBR + P_switch + P_D + P_L;

efficiency = Pout./(Pout+P_losses);

% Efficiency curve
figure;
plot(Pout,efficiency*100);
xlabel('Output power (W)');
ylabel('Efficiency (%)');

% Full load pie chart
X = [P_DBR(1) P_switch(1) P_D(1) P_L(1)];
figure;
p = pie(X);

pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Rectifier losses: ';'MOSFET losses: ';'Diode losses: ';'Inductor losses: '}; 
combinedtxt = strcat(txt,percentValues);

pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
pText(4).String = combinedtxt(4);