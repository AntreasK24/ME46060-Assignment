clc; clear; close all;

%% Objective function
f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;

gradf = @(x) [
    4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7);
    2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)
];

f_grid = @(X,Y) (X.^2 + Y - 11).^2 + (X + Y.^2 - 7).^2;

%% Domain
x = linspace(-6, 6, 400);
y = linspace(-6, 6, 400);
[X, Y] = meshgrid(x, y);
Z = f_grid(X,Y);

%% Initial guess
x0 = [5; 5];

%% Constraints (same structure as SA)
g = {
    @(x) x(1) + x(2) - 1;
    @(x) x(1) - x(2) - 2;   % fixed for feasibility
    @(x) x(1) - 0.2;
    @(x) x(2) + 5;
};

h = {
    
};

%% Jacobians
jaccon = @(x) deal( ...
    [  % Jineq
        1, 1;
        1, -1;
        1, 0;
        0, 1
    ], ...
    [  % Jeq
        1, 0
    ] ...
);

%% Constraint function
confun = @(x) deal( ...
    cellfun(@(gi) gi(x), g)', ...
    cellfun(@(hi) hi(x), h)' ...
);

%% Run SQP
[x_best, cost_best] = SQP(f, gradf, x0, confun, jaccon);

%% =======================
%% PLOTTING
%% =======================

figure;
contour(X, Y, Z, 50);
hold on

%% --- Feasible region (inequalities only) ---
feasible = true(size(X));

for i = 1:length(g)
    Gi = arrayfun(@(x,y) g{i}([x; y]), X, Y);
    feasible = feasible & (Gi <= 0);
end

contourf(X, Y, feasible, 1, ...
    'FaceAlpha', 0.2, ...
    'LineColor', 'none');

%% --- Plot constraints ---
colors = lines(length(g));

labels = {
    '$x_1 + x_2 = 1$'
    '$x_1 - x_2 - 2 = 0$'
    '$x_1 = 0.2$'
    '$x_2 = -5$'
};

for i = 1:length(g)
    Gi = arrayfun(@(x,y) g{i}([x; y]), X, Y);

    contour(X, Y, Gi, [0 0], ...
        'Color', colors(i,:), ...
        'LineWidth', 2);

    text(-5 + i*2, 5 - i*2, labels{i}, ...
        'Interpreter','latex', ...
        'Color', colors(i,:), ...
        'FontSize', 14);
end

%% --- Equality constraint ---
for j = 1:length(h)
    Hj = arrayfun(@(x,y) h{j}([x; y]), X, Y);
    contour(X, Y, Hj, [0 0], 'k--', 'LineWidth', 2);
end

%% --- Solution ---
plot(x_best(1), x_best(2), 'ro', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'r');

%% --- Info box ---
str = sprintf(['SQP solution:\n' ...
               'x_1 = %.3f\n' ...
               'x_2 = %.3f\n' ...
               'f(x) = %.3f'], ...
               x_best(1), x_best(2), cost_best);

annotation('textbox', [0.62, 0.65, 0.3, 0.2], ...
    'String', str, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');

xlabel('$x_1$', 'Interpreter','latex');
ylabel('$x_2$', 'Interpreter','latex');
title('SQP Optimization with Constraints', 'Interpreter','latex');
colorbar

hold off