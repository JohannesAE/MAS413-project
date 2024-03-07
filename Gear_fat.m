clc, clear all, close all
F_t_12 = 1462.99;
F_t_34 = 7014.87;
F_r_12 = 551.27;
F_r_34 = 2553.20;

r1 = 0.04502;
r2 = 0.18221;
r3 = 0.038;
r4 = 0.170;

m = 4;
m_t = 4.14;

% gears
n1 = 21;
n2 = 85;
n3 = 19;
n4 = 85;

% Gear ratios
i1 = n2/n1;
i2 = n4/n3;

i_tot = i1 * i2

%Gear 1 
J_1 = 0.34;
J_2 = 0.41;
J_3 = 1;
J_4 = 1;

Qv = 7;

V_t1 = ((1450*2*pi)/60)*r1
V_t2 = (((1450/i1)*2*pi)/60)*r2
V_t3 = (((1450/i1)*2*pi)/60)*r3
V_t4 = (((1450/i_tot)*2*pi)/60)*r4

K_v1 = 0.72;
K_v2 = 0.72;
K_v3 = 0.875;
K_v4 = 0.875;



K_r = 1.25;
F = 40; %mm



K_a =1.25;
K_m = 1.6;
K_s = 1;
K_B = 1;
K_I = 1;

%Expected life
format bank

life = 10000*1450*60; %10000 hours
K_L = 0.95;
K_T = 1;

%Bending on tooth gear 1
Sig_b1 = (F_t_12/(F*m_t*J_1)) * ((K_a*K_m)/K_v1) *K_s *K_B *K_I;
%Bending on tooth gear 2
Sig_b2 = (F_t_12/(F*m_t*J_2)) * ((K_a*K_m)/K_v2) *K_s *K_B *K_I;
%Bending on tooth gear 3
Sig_b3 = (F_t_34/(F*m_t*J_3)) * ((K_a*K_m)/K_v3) *K_s *K_B *K_I;
%Bending on tooth gear 4
Sig_b4 = (F_t_34/(F*m_t*J_4)) * ((K_a*K_m)/K_v4) *K_s *K_B *K_I;

% Corrected Stress

S_f_uc =  290; %MPA uncorrected fatigue strength for Steel Nitralloy 135M

S_f = (K_L/(K_T*K_r))*S_f_uc;

%Safetyfactors
S_f1 = S_f/Sig_b1
S_f2 = S_f/Sig_b2
S_f3 = S_f/Sig_b3
S_f4 = S_f/Sig_b4