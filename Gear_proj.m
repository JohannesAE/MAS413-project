clc, clear all, close all
%% Parameteres

% gears teeth
n1 = 21;                % Helix
n2 = 85;                % Helix
n3 = 19;                % Spur
n4 = 85;                % Spur

% Gear data
alpha = 20;             % Pressure angle
beta = 15;              % Helix angle
m1 = 4;                 % Module for 1. step
m2 = 4;                 % Module for 2. step
a12 = m1;               % addendum top height
a34 = m2;               % addendum top height
b1 = m1*1.25;           % dedendum

%% Size calculations

% Calculating pitch diameter max 500mm
d1 = (m1*n1)/cosd(beta);                % diameter gear 1
d2 = (m1*n2)/cosd(beta);                % diameter gear 2
d3 = m1*n3;                             % diameter gear 3
d4 = m1*n4;                             % diameter gear 4

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

C1 = (d1+d2)/2;                         % Center distance helical gear ?
C2 = (d3+d4)/2;                         % center distance spur gear

%% Computations

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