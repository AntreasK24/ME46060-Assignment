clc;clear;

Atotal = 200;
Amax = 40;
dT = 18;
Isol = 300;
Uwindow = 1.2;
H = 2500;
ce =2.30;
ct = 120;
cw = 150;

f = @(x) ce*H*( ...
    ((Atotal - x(2))*dT)/(x(3) + x(1)/x(4)) ...
    + Uwindow*x(2)*dT ...
    - x(5)*x(2)*Isol ) ...
    + ct*x(1)*(Atotal - x(2)) ...
    + cw*x(2);

gradf = @(x) [

    % dC/dt
    -ce*H*dT*(Atotal - x(2)) / (x(4)*(x(3) + x(1)/x(4))^2) ...
    + ct*(Atotal - x(2));

    % dC/dAw
    -ce*H*dT/(x(3) + x(1)/x(4)) ...
    + ce*H*Uwindow*dT ...
    - ce*H*x(5)*Isol ...
    - ct*x(1) ...
    + cw;

    % dC/dRbase
    -ce*H*dT*(Atotal - x(2)) / (x(3) + x(1)/x(4))^2;

    % dC/dk
    ce*H*dT*x(1)*(Atotal - x(2)) / ...
    (x(4)^2 * (x(3) + x(1)/x(4))^2);

    % dC/dg
    -ce*H*x(2)*Isol
];

g = {
    @(x) x(1) - 0.15;   % t <= 0.15
    @(x) 0.02 - x(1);   % t >= 0.02 

    @(x) x(2) - 40;     % Aw <= 40
    @(x) 5 - x(2);      % Aw >= 5

    @(x) x(3) - 0.8;    % Rbase <= 0.8
    @(x) 0.4 - x(3);    % Rbase >= 0.4

    @(x) x(4) - 0.04;   % k <= 0.04
    @(x) 0.022 - x(4);   % k >= 0.022

    @(x) x(5) - 0.60;   % g <= 0.60
    @(x) 0.15 - x(5);    % g >= 0.15
};

h = {

};


confun = @(x) deal( ...
    reshape(cellfun(@(gi) gi(x), g), [], 1), ...
    [] ...
);


jaccon = @(x) deal( ...
    [ 1  0  0  0  0;
     -1  0  0  0  0;

      0  1  0  0  0;
      0 -1  0  0  0;

      0  0  1  0  0;
      0  0 -1  0  0;

      0  0  0  1  0;
      0  0  0 -1  0;

      0  0  0  0  1;
      0  0  0  0 -1 ], ...
    [] ...
);

x0 = [0.08; 20; 0.6; 0.03; 0.4];

[x_sa, cost_sa, ~] = simulated_annealing( ...
    f, x0, 100, 1e-3, 0.9, 3000, 0.05, g, h,1e12,true);

[Jineq, ~] = jaccon(x_sa);
[cineq, ~] = confun(x_sa);

disp('Size Jineq:'), disp(size(Jineq))
disp('Size cineq:'), disp(size(cineq))

disp('--- SA Result ---')
disp(x_sa)

[x_sqp, cost_sqp] = SQP(f, gradf, x_sa, confun, jaccon,false,1e8);

disp('--- SQP Result ---')
disp(x_sqp)
disp(cost_sqp)


%Check if KKT conditions are satisified
g = gradf(x_sqp);

mu = zeros(10,1);
mu(1) = -g(1);   % t upper
mu(3) = -g(2);   % Aw upper
mu(5) = -g(3);   % R upper
mu(8) =  g(4);   % k lower
mu(9) = -g(5);   % g upper

[cineq, ceq] = confun(x_sqp);
[Jineq, Jeq] = jaccon(x_sqp);


stationarity = g + Jineq' * mu;
dual_feas    = min(mu);
comp         = mu .* cineq;

disp('stationarity ='), disp(stationarity)
disp('max(cineq,0) ='), disp(max(cineq,0))
disp('dual min     ='), disp(dual_feas)
disp('comp         ='), disp(comp)

% Sensitivity analysis at optimum
grad_opt = gradf(x_sqp);

disp('--- Sensitivity Analysis (gradient at optimum) ---')
vars = {'t','Aw','Rbase','k','g'};

for i = 1:5
    fprintf('dC/d%s = %.4e\n', vars{i}, grad_opt(i));
end

disp('--- Finite-difference sensitivity at optimum ---')

x_ref = x_sqp;
f_ref = f(x_ref);

vars = {'t','Aw','Rbase','k','g'};
rel_step = 1e-4;

for i = 1:5
    x_pert = x_ref;
    h = rel_step * max(abs(x_ref(i)), 1);
    x_pert(i) = x_pert(i) + h;

    sens_fd = (f(x_pert) - f_ref) / h;

    fprintf('FD sensitivity dC/d%s ≈ %.4e\n', vars{i}, sens_fd);
end

disp('--- Normalized sensitivity at optimum ---')

grad_opt = gradf(x_sqp);
f_opt = f(x_sqp);
vars = {'t','Aw','Rbase','k','g'};

for i = 1:5
    S = (x_sqp(i) / f_opt) * grad_opt(i);
    fprintf('Normalized sensitivity for %s = %.4e\n', vars{i}, S);
end

% Bar plot of normalized sensitivities
figure;

vars = {'t','Aw','Rbase','k','g'};
norm_sens = zeros(5,1);

for i = 1:5
    norm_sens(i) = (x_sqp(i) / f(x_sqp)) * grad_opt(i);
end

bar(norm_sens)

set(gca, 'XTickLabel', vars)
xlabel('Design Variables')
ylabel('Normalized Sensitivity')
title('Sensitivity Analysis (Normalized)')
grid on