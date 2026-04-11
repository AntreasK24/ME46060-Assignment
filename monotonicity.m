clear; clc;

% Constants
Atotal  = 200;
dT      = 18;
Uwindow = 1.2;
Isol    = 300;
H       = 2500;
ce      = 2.30;
ct      = 120;
cw      = 150;

% Fix other variables (choose representative values)
Rbase = 0.8;
k     = 0.05;
g     = 0.5;

% Grid
t_vals  = linspace(0.02, 0.30, 50);
Aw_vals = linspace(5, 40, 50);

[T, AW] = meshgrid(t_vals, Aw_vals);

% Compute C
Q = ((Atotal - AW)*dT)./(Rbase + T./k) + Uwindow*AW*dT;
Qsolar = g*AW*Isol;
Qnet = Q - Qsolar;
C = ce*H*Qnet + ct*T.*(Atotal - AW) + cw*AW;

%Surface plot
figure;
surf(T, AW, C)
xlabel('t (insulation thickness)')
ylabel('Aw (window area)')
zlabel('Cost C')
title('Cost function C(t, Aw)')
shading interp
colorbar

%Contour plot 
figure;
contourf(T, AW, C, 20)
xlabel('t')
ylabel('Aw')
title('Contour of C(t, Aw)')
colorbar