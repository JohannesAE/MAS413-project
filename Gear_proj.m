clc, clear all, close all
%% General parameters
nin = 1450;     % input speed
nout = 80.1;    % Output speed
Pin = 8000;     % [kw]
ncoef = 1.29;

% Operating hours
oHours = 10;
days = 365;
year = 4;


%% Parameteres gear

% gears teeth
n1 = 21;                % Helix
n2 = 85;                % Helix
n3 = 19;                % Spur
n4 = 85;                % Spur

% Gear data
alpha = 20;             % Pressure angle
beta = 15;              % Helix angle
m1 = 4/cosd(beta);      % Module for 1. step
m2 = 4;                 % Module for 2. step
a12 = m1;               % addendum top height
a34 = m2;               % addendum top height
b1 = m1*1.25;           % dedendum

%% parameters Shaft/pressfit
d_nom = 25; %mm nominal diameter shaft
%Tolerance params 
L_gear = 40;%mm width of gear sleave
r_nom = 25/2; %Nominal/common radius of fitings
E = 210*1e3;% Gpa to Mpa youngsmodulus plain carbon steel st355
gamma = 0.3; %possions ratio for steel
my_s = 0.2; %static coeff of friction
sigma_y = 355; %Mpa st355 yield strength

% Velocities
nMiddle = nin/(n2/n1);

%% Parameteres bearings
X0 = 1;                     % Radial Factor
Y0 = 1;                     % Axial Factor
X = 1;                      % Dynamic radial Factor
Y = 1;                      % Dynamic Axial Factor
C = 1;                      % Dynamic load ratio[N]
p = 3;                      % Exponent depending on the type of bearing ball=3, roller = 10/3
L10h = oHours*days*year ;   % Operating hours 

%% Parameteres shaft
d_shaft = 0.02;         % Diameter of the shaft

%% Size calculations

% Calculating pitch diameter max 500mm
d1 = (m1*n1)/cosd(beta)                % diameter gear 1
d2 = (m1*n2)/cosd(beta);                % diameter gear 2
d3 = m1*n3;                             % diameter gear 3
d4 = m1*n4;                             % diameter gear 4

% Calculating radius
r1 = d1/2;
r2 = d2/2;
r3 = d3/2;
r4 = d4/3;

% Base circle diameter
db1 = d1-2*b1;                          % Helical 1 ?
db2 = d2-2*b1;                          % Helical 2 ?
db3 = d3*cosd(alpha)                    % Spur gear 3
db4 = d4*cosd(alpha)                    % Spur gear 4 

% Addendum diameter
da1 = d1 + 2*(a12);                     % Helical 1 ?
da2 = d2 + 2*(a12);                     % Helical 2 ?
da3 = d3 + 2*(a34)                      % Spur gear 3
da4 = d4 + 2*(a34)                      % Spur gear 4

C1 = (d1+d2)/2;                         % Center distance helical gear 
C2 = (d3+d4)/2;                         % center distance spur gear

%% Computations gear

% Gear ratios
i1 = n2/n1;
i2 = n4/n3;

i_tot = i1 * i2;

% Calculating contact ratio [CR]
ra1 = da1/2;                            % Addendum radius helical gear 1
ra2 = da2/2;                            % Addendum radius helical gear 2
ra3 = da3/2;                            % Addendum radius spur gear 3
ra4 = da4/2;                            % Addendum radius spur gear 4

rb1 = db1/2;                            % Base circle radius helical gear 1                                                           
rb2 = db2/2;                            % Base circle radius helical gear 2
rb3 = db3/2;                            % Base circle radius spur gear 3
rb4 = db4/2;                            % Base circle radius spur gear 4

pb1 = pi*m1*cosd(alpha);                 % Base pitch ? m1 = mt?
pb2 = pi*m2*cosd(alpha);                % Base pitch

alphawt = acosd((db1+db2)/(C1*2))

num1 = sqrt(ra1^2-rb1^2) + sqrt(ra2^2-rb2^2) - C1*sind(alphawt)
den1 = pb1
CR1 = num1/den1

num2 = sqrt(ra3^2-rb3^2) + sqrt(ra4^2-rb4^2) - C2*sind(alpha)
den2 = pb2
CR2 = num2/den2

%% force Computations 

% calculating load to rear gears       
T_in = (Pin*ncoef)/(2*pi*1450);
T_shaft2 = T_in * (n2/n1);
T_out = T_in * (i_tot);

% Helical gears n1 and n2
F_t_12 = (T_in*10^3)/(r1);
F_r_12 = (F_t_12*tand(alpha))/cosd(beta);
F_a_12 = F_t_12*tand(beta);

% Spur gears n3 and n4
F_t_34 = (T_shaft2*10^3 )/ r3;
F_r_34 = F_t_34 * tand(20);

%% Press fit calculations

%Press fit force
% my_d = 0.90; %Dynamic friction coeff. verdier hentet fra power point
% my_s = 0.10; %Static friction coeff
% p = 0; %Surface pressure
% d = 0; %Common diameter
% 
% F = my_d*p*d*L_gear*pi; %Press force
% 
% T_fit = my_s*p*pi*d*d/2*L_gear; %Torque transmitted from hub to shaft.
% 
% p_min1 = T_in/(my_s*pi*d*d/2*L_gear); %calculating min surface pressure
% 
% p_min2 = T_shaft2/(my_s*pi*d*d/2*L_gear);
% 
% p_min3 = T_out / (my_s*pi*d*d/2*L_gear);

%tolerance
%Choosing a tolerance grade H7/s6 as these tolerances give a snug fit.
%Later calculations will show if it is satisfactory in transmitting torque

%With nominal diameter=25mm and hole tolerance H7/s6 we get an upperlimit
%of 25+0.025mm and lower limit 25+0.0mm
%for the shaft tolerance we get upper limit 25+0.048mm and lower limit
%25+0.035mm
%same tolerance choosen for each shaft

%interference delta
delta_max = (d_nom+0.048)-(d_nom+0);
delta_min = (d_nom+0.035)-(d_nom+0.021);

%Shaft 1%

r_o1 = r1;
r_i1 = 0;

%Surface pressure
sp1_min = (0.5*delta_min)/( (r_nom/E) * ((r_o1^2+r_nom^2)/(r_o1^2-r_i1^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i1^2)/(r_nom^2-r_i1^2))-gamma));
%
sp1_max = (0.5*delta_max)/( (r_nom/E) * ((r_o1^2+r_nom^2)/(r_o1^2-r_i1^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i1^2)/(r_nom^2-r_i1^2))-gamma));

%Transmitable torque
Tt1_min = 2*pi*r_nom^2*sp1_min*L_gear*1e-3;

Nslip1 = Tt1_min/T_in %safety factor against slip

%Stresses/safety factor
sigma_rs1 = -sp1_max;

sigma_th1 = sp1_max*(r_o1^2+r_nom^2)/(r_o1^2-r_nom^2);

sf1_r = sigma_y/sigma_rs1
sf1_t = sigma_y/sigma_th1

%Shaft 2%

r_o2 = r2;
r_i2 = 0;

%Surface pressure
sp2_min = (0.5*delta_min)/( (r_nom/E) * ((r_o2^2+r_nom^2)/(r_o2^2-r_i2^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i2^2)/(r_nom^2-r_i2^2))-gamma));
%
sp2_max = (0.5*delta_max)/( (r_nom/E) * ((r_o2^2+r_nom^2)/(r_o2^2-r_i2^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i2^2)/(r_nom^2-r_i2^2))-gamma));

%Transmitable torque
Tt2_min = 2*pi*r_nom^2*sp2_min*L_gear*1e-3;

Nslip2 = Tt2_min/T_shaft2 %safety factor against slip

%Stresses/safety factor
sigma_rs2 = -sp2_max;

sigma_th2 = sp2_max*(r_o2^2+r_nom^2)/(r_o2^2-r_nom^2);

sf2_r = sigma_y/sigma_rs2
sf2_t = sigma_y/sigma_th2

%Shaft 3%

r_o3 = r3;
r_i3 = 0;

%Surface pressure
sp3_min = (0.5*delta_min)/( (r_nom/E) * ((r_o3^2+r_nom^2)/(r_o3^2-r_i3^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i3^2)/(r_nom^2-r_i3^2))-gamma));
%
sp3_max = (0.5*delta_max)/( (r_nom/E) * ((r_o3^2+r_nom^2)/(r_o3^2-r_i3^2)+gamma) + (r_nom/E) * (((r_nom^2+r_i3^2)/(r_nom^2-r_i3^2))-gamma));

%Transmitable torque
Tt3_min = 2*pi*r_nom^2*sp3_min*L_gear*1e-3;

Nslip3 = Tt3_min/T_out %safety factor against slip

%Stresses/safety factor
sigma_rs3 = -sp3_max;

sigma_th3 = sp3_max*(r_o3^2+r_nom^2)/(r_o3^2-r_nom^2);

sf3_r = sigma_y/sigma_rs3
sf3_t = sigma_y/sigma_th3

%% Bearing computations

% The function of the bearing
% The service life of the bearing
% Calculate the force that occur in the bearing and the direction
% Select the correct bearing type based on the forces acting and the direction of these
% Choose the right bearing size based on the size of the loads and dimensions, e.g. on the shafts

% Selection of bearing type
    % Foreløpig: Tapered roller bearings
    %            Cylnindrical: Radial force dominent, and high
    %            Deep grove: radial, moderate axial, high speed, low noise

% Calculating static bearing load rpm < 10
Po = X0*Fr_bearing + Y0*Fa_bearing ;     % Static bearing load
So = C0/Po;                              % S0 = 0,4 (low) S0 = 4 (high)

% Calculating dynamix bearing load
P_bearing1 = X*Fr_bearing + Y*Fa_bearing;
P_bearing2 = 0;
P_bearing3 = 0;

% Lifetime fo a bearing [24mnd FV]
L10_1 = (L10h*60*nin)/10^6;             % Bearing 1/2
L10_2 = (L10h*60*nMiddle)/10^6          % Bearing 3/4
L10_3 = (L10h*60*nout)/10^6;            % Bearing 5/6                        
C1_bearing = P_bearing1*L10_1^(1/3);
C2_bearing = P_bearing2*L10^(1/3);
C3_bearing = P_bearing3*L10^(1/3);


%%
% Dimensjon shaft
% Hvor skal vi ta opp kreftene
% Hvordan ser akslingene ut (lengde)
% Velge lager
% Formening om hvordan lager festes i boksen
% Tegning diagram
% Ift pitch diameter
