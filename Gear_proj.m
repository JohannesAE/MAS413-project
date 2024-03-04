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

%Shaft
ds1 = 25;%mm diameter shafts

%Tolerance params 
L_gear = 20;%mm width of gear sleave
rn = 25; %Nominal/common radius of fitings
ro = 40; %Outer diameter hub
E = 215;% Mpa youngsmodulus plain carbon steel
gamma = 0.3; %possions ratio for steel
% Velocities
nMiddle = nin/(n2/n1);

%% Parameteres bearings

askf = 0.62;                % This number corresponds to 95% reliability

X1 = 1;                      % Dynamic radial Factor (SKF)
X2 = 1;
X3 = 1;
X4 = 1;
X5 = 1;
X6 = 1;

Y1 = 1;                      % Dynamic Axial Factor(SKF)
Y2 = 1; 
Y3 = 1; 
Y4 = 1;
Y5 = 1;
Y6 = 1;

%Speed for the following bearings
n_bearing12 = nin;          % [rpm]
n_bearing34 = nMiddle;      % [rpm]
n_bearing45 = nout;         % [rpm]

pR = 10/3;                  % Exponent depending on the type of bearing ball=3, roller = 10/3
pB = 3;                     % Exponent bearing

L10h = oHours*days*year ;   % Operating hours 

% Bearings
Bearing1 = 1;               % Ball bearing 
Bearing2 = 2;               % Ball bearing       
Bearing3 = 3;               % Ball bearing
Bearing4 = 4;               % Ball bearing 
Bearing5 = 5;               % Radial load, roller bearing, cap.14.7.3
Bearing6 = 6;               % Radial load, roller bearing

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
%Choosing a tolerance grade H7, h7 and IT7 as these tolerances give a snug fit.
%Later calculations will show if it is satisfactory in transmitting torque
%with shaft diameter 25mm this gives us from table provided in mas412
%25mm+0.025mm upper limit and lower limit 25mm+0.00mm
%and it grade gives a tolerance of 

F_t_34 = (T_shaft2*10^3)/r3;
F_r_34 = F_t_34 * tand(alpha)

%% Bearing computations

% The function of the bearing
% The service life of the bearing
% Calculate the force that occur in the bearing and the direction
% Select the correct bearing type based on the forces acting and the direction of these
% Choose the right bearing size based on the size of the loads and dimensions, e.g. on the shafts

% Choosen Deep grove for 1-4, and roller bearing 5/6

% Forces to the bearings: 
Fr_bearing1 = 1; 
Fa_bearing1 = 1; 
Fr_bearing2 = 1;
Fa_bearing2 = 1; 
Fr_bearing3 = 1; 
Fa_bearing3 = 1;
Fr_bearing4 = 1; 
Fa_bearing4 = 1; 
Fr_bearing5 = 1; 
Fa_bearing5 = 1; 
Fr_bearing6 = 1; 
Fa_bearing6 = 1; 

% Life requirements
L5_1 = (L10h*60*nin)/10^6;             % Bearing 1/2
L5_2 = (L10h*60*nMiddle)/10^6;         % Bearing 3/4
L5_3 = (L10h*60*nout)/10^6;            % Bearing 5/6 

% Compensating for realibility
L10_1 = L5_1/(askf);
L10_2 = L5_2(askf);
L10_3 = L5_3/(askf);

% Calculating dynamix bearing load
P_bearing1 = X1*Fr_bearing1 + Y1*Fa_bearing1;
P_bearing2 = X2*Fr_bearing2 + Y2*Fa_bearing2;
P_bearing3 = X3*Fr_bearing3 + Y3*Fa_bearing3;
P_bearing4 = X4*Fr_bearing4 + Y4*Fa_bearing4;
P_bearing5 = X5*Fr_bearing5 + Y5*Fa_bearing5;
P_bearing6 = X6*Fr_bearing6 + Y6*Fa_bearing6;

% Calculate bearing life
Creg1 = P_bearing1*(L10_1)^(1/pB);       % Bearing 1
Creg2 = P_bearing2*(L10_1)^(1/pB);       % Bearing 2
Creg3 = P_bearing3*(L10_2)^(1/pB);       % Bearing 3
Creg4 = P_bearing4*(L10_2)^(1/pB);       % Bearing 4
Creg5 = P_bearing5*(L10_3)^(1/pR);       % Bearing 5
Creg6 = P_bearing6*(L10_3)^(1/pR);       % Bearing 6

% Now use the skf to find e bearing that satifies Creg
            
%%
% Dimensjon shaft
% Hvor skal vi ta opp kreftene
% Hvordan ser akslingene ut (lengde)
% Velge lager
% Formening om hvordan lager festes i boksen
% Tegning diagram
% Ift pitch diameter
