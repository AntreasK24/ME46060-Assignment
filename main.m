clc; clear; close all;

% Himmelblau's function (anonymous function)
f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
f_grid = @(X,Y) (X.^2 + Y - 11).^2 + (X + Y.^2 - 7).^2;


% Domain
x = linspace(-6, 6, 400);
y = linspace(-6, 6, 400);
[X, Y] = meshgrid(x, y);

x_vec = [X,Y];

% Evaluate function
Z = f_grid(X,Y);

x0 = [5; 5];      % initial point
T0 = 100;         % initial temperature
Tmin = 1e-3;      % stopping temperature
alpha = 0.3;      % cooling rate
maxIter = 1000;
epsilon = 1;

[x_best,cost_best,T] = simulated_annealing(f, x0, T0, Tmin, alpha, maxIter, epsilon);

disp(x_best)

disp(cost_best)

fprintf("Final Temp: %3.2f",T);


% 3D surface plot

z_best = f(x_best);


figure;
surf(X, Y, Z, 'EdgeColor', 'none');
colormap turbo
colorbar
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Himmelblau''s Function');

% Improve visualization
view(45, 45);

hold on
z_best = f(x_best);
plot3(x_best(1), x_best(2), z_best, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off


figure;
contour(X, Y, Z, 50);
hold on
plot(x_best(1), x_best(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('y');
title('Himmelblau''s Function – Contours');
colorbar
hold off
