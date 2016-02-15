%This matlab code takes in the dynamics of the MIP inorder to solve for the
%transfer function to stabilize the mip. The transfer function will then be
%taken used in the code used in the beagle bone black.
close all
clear

Vb = 7.4; %V
wf = 160; %rad/s at Vb
s_bar = .003; %Nm at Vb
Gr = 35.57; %
Im =  3.6*10^-8; % Kg * m^2
Rw = 34e-3; % mm
Mw = 27e-3; % grams each
Mb = 263e-3; % g
L = 36e-3; %mm
Ib = .0004; %Kg * m^2
g = 9.8; % m/s^2

k = .04

Iw = 2*(((Mw*Rw^2)/2) + (Gr^2) * Im)

a = Iw + (Mb + Mw)*Rw^2;
b = Mb*Rw*L;
c = Ib + Mb*L^2;
d = Mb * g * L;

s = tf('s');
%my constants
c1 = (a*c -b^2) / (-a-b);
c2 = (a*d) / (a+b);
b1 = 2*Gr*s_bar;
b2 = 2*(Gr^2)*k;
c3 = (a*c -b^2) / (c+b);
c4 = (b*d) / (c+b);
A1 = (c3*s^2 + b2*s);
A2 = -b2*s;
A4 = (-c4+b2);
A3 = c1*s - c2+b2;

%equation for tf of G1(s)
G1 = (A1*b1 + A2*b1) / (A1*c1*s^2 + A1*c2 - A1*b2*s + A2*c4 - A2*b2*s);
G1 = minreal(G1)
num1 = -706.8*s
den1 = s^3 + 7.923e05*s^2 - 180.4*s - 3.991e07

%equation for tf of G2(s)
G2 = (-A4*b1 + A3*b1) / (A3*c3*s^2 + A3*b2*s - A4*b2*s);
G2 = minreal(G2)
num2 = -(A4*b1 - A3*b1);
num2 = minreal(num2);
den2 = (A3*c3*s^2 - A3*b2*s - A4*b2*s);
den2 = minreal(den2);

% Graphs of Impulse responce to each
figure
impulse(G1)
title('Impulse responce to G1(s)')
axis([0 1.2 -1 .4])

figure
impulse(G2)
title('Impulse responce to G2(s)')

[A,B,C,D] = tf2ss([-706.8 0],[1 7.923e5 -180.4 3.991e7]);
eig(A)
[z,p,k] = zpkdata(G1)
%Bode and root locus to start stabilizing G1(s)
figure
rlocus(1/G1)
figure
bode(G1), grid
title('bode G1')

%design D1(s)
D1 = (s + .1) / (s+10)
K = .042
figure
rlocus(D1*G1)
title('D1*G1')

figure
rlocus((D1*G1) / (1 + D1*G1))
title('closed loop (D1*G1) / (1 + D1*G1')

figure
bode((D1*G1) / (1 + D1*G1)), grid

%put into discrete time using tustins approx
%G1_d = c2d(G1,.005,'tustin')
%G2_d = c2d(G2,.005,'tustin')

G11 = tf([-669.4 0] , [1 31.26 -175.6 -1549]);% actual plant equation
D11 = tf([1 5],[1 0]); %filter should work
