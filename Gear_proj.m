 clc, clear all, close all

% gears teeth
n1 = 21;
n2 = 85;
n3 = 19;
n4 = 85;

% Gear ratios
i1 = n2/n1;
i2 = n4/n3;

i_tot = i1 * i2

% gear data
m1 = 4;                 %Module for 1. step
a1 = m1;                %addendum
b1 = m1*1.25;           %dedendum

% Calculating pitch diameter max 500mm
d1 = m1*n1;              %diameter gear 1
d2 = m1*n2;              %diameter gear 2
d3 = m1*n3;              %diameter gear 3
d4 = m1*n4;              %diameter gear 4

% Calculating contact ratio [CR]
CR = 

r1 = d1/2;
r2 = d2/2;
C = r1+r2;              %center distance

bcr_1 = r1*cosd(20)     %base cicrle radius
bcr_2 = r2*cosd(20)

pc1 = (pi*d1)/n1        %pitch circle g1
pc2 = (pi*d2)/n2

pb1 = pc1*cosd(20)
pb2 = pc2*cosd(20)

m2 = 6;
