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



% Inequality constraints g(x) <= 0
g = {
    @(x) x(1) + x(2) - 1;        % x1 + x2 <= 3
    @(x) x(1) + x(2) + 1;

};

% Equality constraints h(x) = 0 
h = {
    % @(x) x(1)^2 + x(2)^2 - 5;
};


[x_best,cost_best,T] = simulated_annealing(f, x0, T0, Tmin, alpha, maxIter, epsilon, g, h);

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

figure;
contour(X, Y, Z, 50);
hold on


plot(x_best(1), x_best(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');


G_total_violation = zeros(size(X)); 

for i = 1:length(g)

    Gi = arrayfun(@(x,y) g{i}([x; y]), X, Y);


    contour(X, Y, Gi, [0 0], 'LineWidth', 2);


    G_total_violation = G_total_violation + (Gi > 0);
end


contourf(X, Y, G_total_violation > 0, 1, ...
    'FaceAlpha', 0.2, 'LineColor', 'none');

xlabel('x');
ylabel('y');
title('Himmelblau''s Function with Constraints');
colorbar
hold off