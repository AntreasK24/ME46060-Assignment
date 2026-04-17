clc;clear;

%Constants

Rbase = 0.6;
k = 0.03;
Atotal = 200;
Amax = 40;
dT = 18;
g_sol = 0.4;
Isol = 300;
Uwindow = 1.2;
H = 2500;
ce =2.30;
ct = 120;
cw = 150;

% Simplified Problem

% Objective
f = @(x) ce*H*( ...
    ((Atotal - x(2))*dT)/(Rbase + x(1)/k) ...
    + Uwindow*x(2)*dT ...  
    - g_sol*x(2)*Isol ) ...
    + ct*x(1)*(Atotal - x(2)) ...
    + cw*x(2);

gradf = @(x) [
    -ce*H*dT*(Atotal - x(2)) / (k*(Rbase + x(1)/k)^2) ...
    + ct*(Atotal - x(2));

    -ce*H*dT/(Rbase + x(1)/k) ...
    + ce*H*Uwindow*dT ...
    - ce*H*g_sol*Isol ...
    - ct*x(1) ...
    + cw
];

g = {
    @(x) x(1) - 0.15;   % t <= 0.15
    @(x) 0.02 - x(1);   % t >= 0.02
    @(x) x(2) - 40;     % Aw <= 40
    @(x) 5 - x(2);      % Aw >= 5
};

h = {};



confun = @(x) deal( ...
    cell2mat(cellfun(@(gi) gi(x), g, 'UniformOutput', false)), ...
    [] ...
);


jaccon = @(x) deal( ...
    [ 1  0;
     -1  0;
      0  1;
      0 -1], ...
    [] ...
);



x0 = [0.1; 20];

[x_sa, cost_sa, ~] = simulated_annealing( ...
    f, x0, 100, 1e-3, 0.9, 3000, 0.001, g, h,8e10,true);

[Jineq, ~] = jaccon(x_sa);
[cineq, ~] = confun(x_sa);

disp('Size of Jineq:'); disp(size(Jineq))
disp('Size of cineq:'); disp(size(cineq))

disp('--- SA Result ---')
disp(x_sa)

x0 = [0.1; 20];
[x_sqp, cost_sqp] = SQP(f, gradf, x_sa, confun, jaccon,false,1e8);


disp('--- SQP Result ---')
disp(x_sqp)


% Expanded grid
t = linspace(0, 0.15, 300);
Aw = linspace(5, 40, 300);
[T, AW] = meshgrid(t, Aw);

% Cost function
Z = ce*H*( ...
    ((Atotal - AW)*dT)./(Rbase + T./k) ...
    + Uwindow*AW*dT ...
    - g_sol*AW*Isol ) ...
    + ct*T.*(Atotal - AW) ...
    + cw*AW;

%Main plot


figure;

fontSize = 20;


contour(T, AW, Z, 50);
hold on

% Constraint boundaries
h1 = plot([0.02 0.15],[5 5],'k','LineWidth',2);
h2 = plot([0.02 0.15],[40 40],'k','LineWidth',2);
h3 = plot([0.02 0.02],[5 40],'k','LineWidth',2);
h4 = plot([0.15 0.15],[5 40],'k','LineWidth',2);

% Plot solutions
h_sa  = plot(x_sa(1), x_sa(2), 'bo','MarkerFaceColor','b','MarkerSize',8);
h_sqp = plot(x_sqp(1), x_sqp(2), 'ro','MarkerFaceColor','r','MarkerSize',10);

xlabel('t','FontSize',fontSize)
ylabel('A_w','FontSize',fontSize)
title('Contour plot with constraints','FontSize',fontSize)
colorbar
grid on

legend([h1 h_sa h_sqp], ...
    {'Constraints','SA','SQP'}, ...
    'Location','northwest','FontSize',fontSize)

% %Zoomed inset
% axes('Position',[0.25 0.2 0.3 0.3]) % position of inset
% box on

% contour(T, AW, Z, 30);
% hold on

% % Replot constraints
% plot([0.02 0.30],[5 5],'k','LineWidth',1.5)
% plot([0.02 0.30],[40 40],'k','LineWidth',1.5)
% plot([0.02 0.02],[5 40],'k','LineWidth',1.5)
% plot([0.30 0.30],[5 40],'k','LineWidth',1.5)

% % Replot points
% plot(x_sa(1), x_sa(2), 'bo','MarkerFaceColor','b')
% plot(x_sqp(1), x_sqp(2), 'ro','MarkerFaceColor','r')

% % Zoom window (tight around solution)
% xlim([0.28 0.305])
% ylim([38 41])

% title('Zoom near optimum','FontSize',fontSize)
% grid on

