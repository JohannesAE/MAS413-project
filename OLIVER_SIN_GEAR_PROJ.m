 clc, clear all, close all

% gears
n1 = 21;
n2 = 85;
n3 = 19;
n4 = 85;

% Gear ratios
i1 = n2/n1;
i2 = n4/n3;

i_tot = i1 * i2

for i = 1
% gear data for spur gears
m1 = 4;     %Module for 1. step
a1 = m1;    %addendum
b1 = m1*1.25;   %dedendum

d3 = m1*n3;      %diameter gear 3
r3 = d3/2;
d4 = m1*n4;    %diameter gear 2
r4 = d4/2
C = r3+r4;      %center distance

bcr_3 = r3*cosd(20); %base cicrle radius
bcr_4 = r4*cosd(20);

pc3 = (pi*d3)/n3; %pitch circle g1
pc4 = (pi*d4)/n4;

pb3 = pc3*cosd(20);
pb4 = pc4*cosd(20);


% Calculating the contact ratio

% Length of action "Z"
Z = sqrt((r3+a1)^2 - (r3*cosd(20))^2) + sqrt((r4+a1)^2 - (r4*cosd(20))^2) - C*sind(20);

%diametral pitch "Pd"
pd = n3/(2*r3);

% contact ratio "mp"
m_p = (pd*Z)/(pi*cosd(20)); %1.69 som ebbesen mener er greit for spur gears



end

% Helical %
alpha = 20;
beta =  15;
m_t = m1/cosd(beta);
a_t = m_t;
b_t = m_t*1.25;

% Pitch diameter
d1 = (m_t *n1)/cosd(beta);
d2 = (m_t *n2)/cosd(beta);
r1 = d1/2;
r2 = d2/2;
C12 = r1+r2;

da1 = d1 + 2*m_t;
db1 = d1 - 2*m_t*1.25;
da2 = d2 + 2*m_t;
db2 = d2 - 2*m_t*1.25;
awt = acosd((db1 + db2)/(2*C12))

CR_helix = ( sqrt((da1/2)^2 - (db1/2)^2) + sqrt( (da2/2)^2 - (db2/2)^2) -C12*sind(awt) )/(pi*m_t*cosd(20))    


           % calculating load to rear gears
Pin = 8000;        
T_in = (Pin*1.25*60)/(2*pi*1450);
T_shaft2 = T_in * (n2/n1);
T_out = T_in * (i_tot);
format bank
% Helical gears n1 and n2
F_t_12 = (T_in*10^3)/(r1)
F_r_12 = (F_t_12*tand(alpha))/cosd(beta)
F_a_12 = F_t_12*tand(beta)

% Spur gears n3 and n4
F_t_34 = (T_shaft2*10^3 )/ r3
F_r_34 = F_t_34 * tand(20)


% SHAFT 1
for i = 0
A_x_xz_1 = F_a_12;
B_y_xz = (F_t_12*0.1)/0.3;
A_y_xz = F_t_12-B_y_xz;

A_x_xy_1 = F_a_12
B_y_xy = ((F_r_12*0.1)+(F_a_12*(0.04502)))/0.3;
A_y_xy = F_r_12-B_y_xy;


x1 = 0 : 0.01 : 0.1;
x11 = 0.1 : 0.01 : 0.3;
%Shaft 1 (XZ PLANE)
p111XZ = F_a_12*ones(size(x1));
V111XZ = -A_y_xz*ones(size(x1));
M111XZ = -A_y_xz*x1;

p112XZ = (A_x_xz_1 - F_a_12)*ones(size(x11));
V112XZ = (F_t_12 - A_y_xz)*ones(size(x11));
M112XZ = F_t_12 *(x11-0.1) - A_y_xz*(x11);



%%%%%%%%%%%% XY PLANE %%%%%%%%%%%%%
p111XY = A_x_xy_1*ones(size(x1));
V111XY = A_y_xy*ones(size(x1));
M111XY = A_y_xy*x1;

p112XY =( A_x_xy_1-F_a_12)*ones(size(x11));
V112XY = (-F_r_12 + A_y_xy)*ones(size(x11));
M112XY = A_y_xy*x11 - (F_r_12*(x11-0.1)) + (F_a_12*(0.04502));
end


%%%%%%%%%%%%%%% Shaft 2 %%%%%%%%%%%%%
for i = 0
Cx_xz = F_a_12;
Dz = (F_t_34*0.2 -(F_t_12*0.1))/0.3;
Cz = F_t_12 - F_t_34 + Dz;

%Region 1 xz plane
x21 = 0: 0.01 : 0.1;
Px_xz_21 = Cx_xz * ones(size(x21));
Vx_xz_21 = Cz * ones(size(x21));
Mx_xz_21 = Cz*x21;
%Region 2 xz plane
x22 = 0.1 : 0.01 : 0.2;
Px_xz_22 = (Cx_xz - F_a_12)*ones(size(x22));
Vx_xz_22 = (Cz-F_t_12)*ones(size(x22));
Mx_xz_22 = Cz*x22 - (F_t_12*(x22-0.1));
%Region 3 xz plane
x23 = 0.2 : 0.01 : 0.3;
Px_xz_23 = (Cx_xz - F_a_12)*ones(size(x23));
Vx_xz_23 = (-F_t_12 + F_t_34 + Cz)* ones(size(x23));
Mx_xz_23 = F_t_34*(x23-0.2) - F_t_12*(x23-0.1) + Cz*(x23);


%%% SHAFT 2 XY PLANE %%%
Cx_xy = F_a_12;
Dy = (-(F_a_12*0.182) - (F_r_12*0.1) + (F_r_34 * 0.2))/0.3;
Cy = F_r_12 + F_r_34 -Dy;
%Region 1
Px_xy_21 = -Cx_xy*ones(size(x21));
Vx_xy_21 = -Cy * ones (size(x21));
Mx_xy_21 = -Cy*x21;
%Region 2 
Px_xy_22 = (Cx_xy - F_a_12)*ones(size(x22));
Vx_xy_22 = (F_r_12-Cy)*ones(size(x22));
Mx_xy_22 = -(F_a_12*0.182) + (F_r_12*(x22-0.1)) - (Cy*(x22));
%Region 3
Px_xy_23 = (Cx_xy - F_a_12)*ones(size(x23));
Vx_xy_23  = (F_r_12-F_r_34 - Cy)*ones(size(x23));
Mx_xy_23 = (F_r_12*(x23-0.1)) - (F_a_12*0.182) -(Cy*x23) - (F_r_34*(x23-0.2));
end



%%%%% SHAFT 3 %%%
for i = 0
%XZ plane
Fz = (F_t_34*0.2)/0.3;
Ez = (F_t_34-Fz);

% Region 1
x31 = 0:0.01:0.2;
Px_xz_31 = zeros(size(x31));
Vx_xz_31 = Ez * ones(size(x31));
Mx_xz_31 = Ez*x31;
% Region 2 
x32 = 0.2:0.01:0.3;
Px_xz_32 = zeros(size(x32));
Vx_xz_32 = (Ez-F_t_34)*ones(size(x32));
Mx_xz_32 = Ez*x32 - (F_t_34*(x32-0.2));

% XY PLANE
Fy = (F_r_34*0.2)/0.3;
Ey = F_r_34-Fy;

%Region 1
Px_xy_31 = zeros(size(x31));
Vx_xy_31 = -Ey*ones(size(x31));
Mx_xy_31 =  -Ey*x31;
% Region 2
Px_xy_32 = zeros(size(x32));
Vx_xy_32 = (F_r_34 -Ey)*ones(size(x32));
Mx_xy_32 = (F_r_34*(x32-0.2)) - (Ey*x32);
end


%%%%%%% PLOTS %%%%%%%%%%%%%%

for i = 0;
% SHAFT 1 
for i = 0;
% XZ PLANE
subplot(3,1,1)
plot(x1,p111XZ,'LineStyle','-')
title('XZ-Plane, SHAFT 1')
hold on
plot(x11,p112XZ,'LineStyle','-')
grid on
subplot(3,1,2)
plot(x1,V111XZ)
hold on
plot(x11,V112XZ)
ylim([350 370])
grid on
subplot(3,1,3)
plot(x1, M111XZ)
hold on
plot(x11,M112XZ)
grid on



%%XY PLANE
figure
subplot(3,1,1)
plot(x1,p111XY)
title('XY-Plane, SHAFT 1')
hold on
plot(x11,p112XY,'LineStyle','-')
grid on
subplot(3,1,2)
plot(x1,V111XY)
hold on
plot(x11,V112XY)

grid on
subplot(3,1,3)
plot(x1, M111XY)
hold on
plot(x11,M112XY)
grid on

end
% SHAFT 2 
for i = 0;
% XZ PLANE
figure
subplot(3,1,1)
plot(x21,Px_xz_21)
title('XZ-Plane, SHAFT 2')
hold on
plot(x22,Px_xz_22,'LineStyle','-')
plot(x23,Px_xz_23,'LineStyle','-')
grid on

subplot(3,1,2)
plot(x21,Vx_xz_21)
hold on
plot(x22,Vx_xz_22)
plot(x23,Vx_xz_23)
%ylim([-300 350])
grid on
subplot(3,1,3)
plot(x21, Mx_xz_21)
hold on
plot(x22, Mx_xz_22)
plot(x23, Mx_xz_23)
grid on
%ylim([-25 25])

% XY PLANE
figure
subplot(3,1,1)
plot(x21,Px_xy_21)
title('XY-Plane, SHAFT 2')
hold on
plot(x22,Px_xy_22,'LineStyle','-')
plot(x23,Px_xy_23,'LineStyle','-')
grid on

subplot(3,1,2)
plot(x21,Vx_xy_21)
hold on
plot(x22,Vx_xy_22)
plot(x23,Vx_xy_23)
%ylim([-300 350])
grid on
subplot(3,1,3)
plot(x21, Mx_xy_21)
hold on
plot(x22, Mx_xy_22)
plot(x23, Mx_xy_23)
grid on
%ylim([-25 25])

end
% SHAFT 3
for i = 0;
% XZ PLANE
figure
subplot(3,1,1)
plot(x31,Px_xz_31)
title('XZ-Plane, SHAFT 3')
hold on
plot(x32,Px_xz_32,'LineStyle','-')
grid on

subplot(3,1,2)
plot(x31,Vx_xz_31)
hold on
plot(x32,Vx_xz_32)
%ylim([-300 350])

grid on
subplot(3,1,3)
plot(x31, Mx_xz_31)
hold on
plot(x32, Mx_xz_32)
grid on
%ylim([-25 25])

% XY PLANE
figure
subplot(3,1,1)
plot(x31,Px_xy_31)
title('XY-Plane, SHAFT 3')
hold on
plot(x32,Px_xy_32,'LineStyle','-')
grid on

subplot(3,1,2)
plot(x31,Vx_xy_31)
hold on
plot(x32,Vx_xy_32)
%ylim([-300 350])

grid on
subplot(3,1,3)
plot(x31, Mx_xy_31)
hold on
plot(x32, Mx_xy_32)
grid on
%ylim([-25 25])
end


end

%Torque between gears on shaft 2
T_diff = T_in*i1 - T_in*i_tot
figure 
plot(x21,zeros(size(x21)))
title('Torque shaft 2 Nm')
hold on
grid on
plot(x22, T_diff*ones(size(x22)))
plot(x23,zeros(size(x23)))


% Bearing loads for shaft 1
A_load = sqrt(A_y_xy^2 + A_y_xz^2)
B_load = sqrt(B_y_xy^2 + B_y_xz^2)
C_load = sqrt(Cy^2+Cz^2)
D_load = sqrt(Dy^2 + Dz^2)
E_load = sqrt(Ey^2 + Ez^2)
F_load = sqrt(Fy^2 + Fz^2)

