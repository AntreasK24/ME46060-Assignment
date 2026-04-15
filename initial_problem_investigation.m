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

syms t Aw Rbase k g real
x = [t, Aw, Rbase, k, g];

% Model
Q = ((Atotal - Aw)*dT)/(Rbase + t/k) + Uwindow*Aw*dT;
Qsolar = g*Aw*Isol;
Qnet = Q - Qsolar;
C = ce*H*Qnet + ct*t*(Atotal - Aw) + cw*Aw;


%% Convexity Analysis

% Hessian of cost function
H_C = simplify(hessian(C, x));
Hfun_C = matlabFunction(H_C, 'Vars', {t, Aw, Rbase, k, g});

% Feasible region
t_vals     = linspace(0.02, 0.15, 6);
Aw_vals    = linspace(5, 40, 6);
Rbase_vals = linspace(0.4, 0.8, 5);
k_vals     = linspace(0.022, 0.04, 5);
g_vals     = linspace(0.15, 0.60, 5);

minEig = inf;
worstPoint = nan(1,5);

for t0 = t_vals
    for Aw0 = Aw_vals
        for R0 = Rbase_vals
            for k0 = k_vals
                for g0 = g_vals
                    Hnum = Hfun_C(t0, Aw0, R0, k0, g0);
                    Hnum = (Hnum + Hnum.')/2;
                    e = eig(Hnum);
                    mine = min(e);

                    if mine < minEig
                        minEig = mine;
                        worstPoint = [t0, Aw0, R0, k0, g0];
                    end
                end
            end
        end
    end
end

fprintf('Smallest eigenvalue found: %.12g\n', minEig);
fprintf('Worst point [t, Aw, Rbase, k, g]: [%.6g, %.6g, %.6g, %.6g, %.6g]\n', ...
    worstPoint(1), worstPoint(2), worstPoint(3), worstPoint(4), worstPoint(5));

if minEig < -1e-8
    fprintf('Conclusion: C is NOT convex on the tested region.\n');
elseif minEig <= 1e-8
    fprintf('Conclusion: C appears convex / PSD on the tested region.\n');
else
    fprintf('Conclusion: C appears strongly convex on the tested region.\n');
end

%% Monotonicity 
% Track min/max derivative values for each variable
minGrad = inf(1,5);
maxGrad = -inf(1,5);
minPoint = nan(5,5);   
maxPoint = nan(5,5);   

gradC = simplify(gradient(C, x));
gradFun = matlabFunction(gradC, 'Vars', {t, Aw, Rbase, k, g});


for t0 = t_vals
    for Aw0 = Aw_vals
        for R0 = Rbase_vals
            for k0 = k_vals
                for g0 = g_vals
                    gnum = gradFun(t0, Aw0, R0, k0, g0);
                    gnum = gnum(:).';   % make row vector

                    point = [t0, Aw0, R0, k0, g0];

                    for i = 1:5
                        if gnum(i) < minGrad(i)
                            minGrad(i) = gnum(i);
                            minPoint(i,:) = point;
                        end
                        if gnum(i) > maxGrad(i)
                            maxGrad(i) = gnum(i);
                            maxPoint(i,:) = point;
                        end
                    end
                end
            end
        end
    end
end

vars = {'t','Aw','Rbase','k','g'};
tol = 1e-8;

fprintf('\nMonotonicity results on tested region:\n');
for i = 1:5
    fprintf('\nVariable %s:\n', vars{i});
    fprintf('  min partial derivative = %.12g\n', minGrad(i));
    fprintf('  at point [%.6g, %.6g, %.6g, %.6g, %.6g]\n', minPoint(i,1), minPoint(i,2), minPoint(i,3), minPoint(i,4), minPoint(i,5));
    fprintf('  max partial derivative = %.12g\n', maxGrad(i));
    fprintf('  at point [%.6g, %.6g, %.6g, %.6g, %.6g]\n', maxPoint(i,1), maxPoint(i,2), maxPoint(i,3), maxPoint(i,4), maxPoint(i,5));

    if minGrad(i) >= -tol
        fprintf('  Conclusion: C is monotonically INCREASING in %s on the tested region.\n', vars{i});
    elseif maxGrad(i) <= tol
        fprintf('  Conclusion: C is monotonically DECREASING in %s on the tested region.\n', vars{i});
    else
        fprintf('  Conclusion: C is NOT monotonic in %s on the tested region.\n', vars{i});
    end
end

% Fix other variables 
Rbase = 0.6;
k     = 0.03;
g     = 0.4;

% Grid
t_vals  = linspace(0.02, 0.15, 50);
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


% Partial minimization over Aw
t_grid = linspace(0.02, 0.15, 100);
phi = zeros(size(t_grid));
Aw_star = zeros(size(t_grid));

for i = 1:length(t_grid)
    t_fixed = t_grid(i);

    cost_Aw = @(Aw) ce*H*(((Atotal - Aw)*dT)/(Rbase + t_fixed/k) + Uwindow*Aw*dT - g*Aw*Isol ) + ct*t_fixed*(Atotal - Aw) + cw*Aw;

    Aw_candidates = linspace(5, 40, 200);
    costs = arrayfun(cost_Aw, Aw_candidates);

    [phi(i), idx] = min(costs);
    Aw_star(i) = Aw_candidates(idx);
end

figure;
plot(t_grid, phi, 'LineWidth', 2);
xlabel('t');
ylabel('\phi(t) = min_{A_w} C(t,A_w)');
title('Partial minimization over window area');
grid on;