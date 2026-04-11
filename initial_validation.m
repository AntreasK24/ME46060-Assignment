clc; clear; close all;

f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;

gradf = @(x) [
    4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7);
    2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)
];


g = {
    @(x) x(1) + x(2) - 1;
    @(x) x(1) - x(2) - 2;
    @(x) x(1) - 0.2;
    @(x) x(2) + 3.5;
};

h = {
  
};

confun = @(x) deal( ...
    cellfun(@(gi) gi(x), g)', ...
    cellfun(@(hi) hi(x), h)' ...
);

jaccon = @(x) deal( ...
    [1 1;
     1 -1;
     1 0;
     0 1], ...
    [] ...
);


x = linspace(-6, 6, 400);
y = linspace(-6, 6, 400);
[X, Y] = meshgrid(x, y);

Z = (X.^2 + Y - 11).^2 + (X + Y.^2 - 7).^2;


x0 = [5; 5];

T0 = 100;
Tmin = 1e-3;
alpha = 0.9;
maxIter = 3000;
epsilon = 1;

[x_sa, cost_sa, ~] = simulated_annealing( ...
    f, x0, T0, Tmin, alpha, maxIter, epsilon, g, h,1000,true);

disp('--- SA Result ---')
disp(x_sa)


[x_sqp, cost_sqp] = SQP(f, gradf, x_sa, confun, jaccon,false,1000);

disp('--- SQP Result ---')
disp(x_sqp)
disp(cost_sqp)

figure;
contour(X, Y, Z, 50);
hold on



colors = lines(length(g));

labels = {
    '$x_1 + x_2 = 1$'
    '$x_1 - x_2 - 2 = 0$'
    '$x_1 = 0.2$'
    '$x_2 = -3.5$'
};

for i = 1:length(g)
    Gi = arrayfun(@(x,y) g{i}([x; y]), X, Y);

    contour(X, Y, Gi, [0 0], ...
        'Color', colors(i,:), ...
        'LineWidth', 2);


    text(-5 + i*2, 5 - i*2, labels{i}, ...
        'Interpreter','latex', ...
        'Color', colors(i,:), ...
        'FontSize', 21, ...
        'FontWeight','bold');
end


h_sa = plot(x_sa(1), x_sa(2), 'bo', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'b');

h_sqp = plot(x_sqp(1), x_sqp(2), 'ro', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'r');

legend([h_sa, h_sqp], {'SA', 'SQP'}, 'Location','northeast')

str = sprintf(['SA → SQP result:\n' ...
               'x_1 = %.3f\n' ...
               'x_2 = %.3f\n' ...
               'f(x) = %.3f'], ...
               x_sqp(1), x_sqp(2), f(x_sqp));

annotation('textbox', [0.15, 0.4, 0.3, 0.2], ...
    'String', str, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w');


xlabel('$x_1$', 'Interpreter','latex');
ylabel('$x_2$', 'Interpreter','latex');
title('Hybrid SA $\rightarrow$ SQP Optimization with Constraints', ...
    'Interpreter','latex');

colorbar
hold off