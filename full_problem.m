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
    @(x) x(1) - 0.30;   % t <= 0.30
    @(x) 0.02 - x(1);   % t >= 0.02 

    @(x) x(2) - 40;     % Aw <= 40
    @(x) 5 - x(2);      % Aw >= 5

    @(x) x(3) - 1.5;    % Rbase <= 1.5
    @(x) 0.2 - x(3);    % Rbase >= 0.2

    @(x) x(4) - 0.08;   % k <= 0.08
    @(x) 0.02 - x(4);   % k >= 0.02

    @(x) x(5) - 0.75;   % g <= 0.75
    @(x) 0.2 - x(5);    % g >= 0.2
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

x0 = [0.1; 20; 0.8; 0.05; 0.5];

[x_sa, cost_sa, ~] = simulated_annealing( ...
    f, x0, 100, 1e-3, 0.9, 3000, 0.05, g, h);

[Jineq, ~] = jaccon(x_sa);
[cineq, ~] = confun(x_sa);

disp('Size Jineq:'), disp(size(Jineq))
disp('Size cineq:'), disp(size(cineq))

disp('--- SA Result ---')
disp(x_sa)

[x_sqp, cost_sqp] = SQP(f, gradf, x_sa, confun, jaccon,false);

disp('--- SQP Result ---')
disp(x_sqp)
disp(cost_sqp)