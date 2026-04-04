clear; clc;

% Objective function
f = @(x) (x(1)-1)^2 + (x(2)-2)^2;

% Gradient of objective
gradf = @(x) [2*(x(1)-1);
              2*(x(2)-2)];

% Constraint function
% inequality: x1^2 + x2^2 - 1 <= 0
% no equality constraints
confun = @(x) deal( ...
    x(1)^2 + x(2)^2 - 1, ... % cineq
    [] ...                   % ceq
);

% Jacobian of constraints
jaccon = @(x) deal( ...
    [2*x(1), 2*x(2)], ... % Jineq (row vector)
    [] ...                % Jeq
);

% Initial guess
x0 = [0.5; 0.5];

% Call your SQP solver
[x_opt, cost_opt] = SQP(f, gradf, x0, confun, jaccon);

% Display results
disp('Optimal x:')
disp(x_opt)

disp('Optimal cost:')
disp(cost_opt)