    % Rbase = 0.4;
    % k = 0.035;
    Atotal = 200;
    Amax = 40;
    dT = 18;
    % g = 0.6;
    Isol = 300;
    Uwindow = 1.2;
    H = 2500;
    ce =2.30;
    ct = 120;
    cw = 150;
    gamma = 0.3;
    betta_cool = 5;

    Rwall = @(t,Rbase,k) Rbase + 2*t/k;

    Qwall = @(t,Aw) (Atotal - Aw)*dT./Rwall(t);

    Qwindow = @(Aw) Uwindow*Aw*dT;

    Qsolar = @(Aw,g) g*Isol*Amax*(1 - exp(-Aw./Amax));

    Qcool = @(t,Aw) betta_cool*Aw.^2 .* (1 + 5*t);

    Qnet = @(t,Aw) Qwall(t,Aw) + Qwindow(Aw) + Qsolar(Aw);

    t_star = 0.18;
    Aw_star = 30;

    lambda_t = 8e7;
    lambda_A = 200;

    Ccomfort = @(t,Aw) ...
        lambda_t*(t - t_star).^2 + ...
        lambda_A*(Aw - Aw_star).^2;




    E = @(t,Aw) H * Qnet(t,Aw);

    C = @(t,Aw) ...
        ce*H*Qnet(t,Aw)/1000 ...
        + ce*H*Qcool(t,Aw)/1000 ...
        + cw*Aw ...
        + gamma*Aw.^2 ...
        + 5e5*(Aw/Atotal - 1.5*t).^2 ...
        + Ccomfort(t,Aw);


    %Sanity Checks

    t = linspace(0.01,0.4,100);

    figure
    plot(t,Rwall(t),'LineWidth',2)
    xlabel('Insulation thickness t (m)')
    ylabel('R_{wall} (m^2K/W)')
    title('Wall Thermal Resistance')
    grid on

    t = linspace(0.01,0.4,100);
    Aw = 40;

    figure
    plot(t,Qwall(t,Aw),'LineWidth',2)
    xlabel('Insulation thickness t (m)')
    ylabel('Q_{wall} (W)')
    title('Wall Heat Loss vs Insulation Thickness')
    grid on

    Aw = linspace(0,80,100);

    figure
    plot(Aw,Qwindow(Aw),'LineWidth',2)
    hold on
    plot(Aw,Qsolar(Aw),'LineWidth',2)
    xlabel('Window area A_w (m^2)')
    ylabel('Heat Flow (W)')
    title('Window Heat Loss vs Solar Gain')
    legend('Window Heat Loss','Solar Gain')
    grid on


    t = linspace(0.01,0.4,100);
    Aw = 40;

    figure
    plot(t,Qnet(t,Aw),'LineWidth',2)
    xlabel('Insulation thickness t (m)')
    ylabel('Q_{net} (W)')
    title('Net Heating Demand vs Insulation Thickness')
    grid on

    Aw = linspace(0,80,100);
    t = 0.15;

    figure
    plot(Aw,Qnet(t,Aw),'LineWidth',2)
    xlabel('Window area A_w (m^2)')
    ylabel('Q_{net} (W)')
    title('Net Heating Demand vs Window Area')
    grid on

    t = linspace(0.01,0.4,100);
    Aw = 40;

    figure
    plot(t,C(t,Aw),'LineWidth',2)
    xlabel('Insulation thickness t (m)')
    ylabel('Cost (€)')
    title('Cost vs Insulation Thickness')
    grid on

    [tg,Awg] = meshgrid(linspace(0.01,0.4,200), linspace(0,80,50));

    Cost = C(tg,Awg);

    figure
    surf(tg,Awg,Cost)
    xlabel('Insulation thickness t')
    ylabel('Window area A_w')
    zlabel('Cost (€)')
    title('Cost Surface')

    figure
    contour(tg,Awg,Cost,30)
    xlabel('Insulation thickness t')
    ylabel('Window area A_w')
    title('Cost Contours')
    colorbar

t = linspace(0.01,0.4,200);
Aw_fixed = 40;

figure
plot(t, C(t,Aw_fixed), 'LineWidth',2)
xlabel('Insulation thickness t (m)')
ylabel('Cost C')
title('Cost vs Insulation Thickness (A_w = 40 m^2)')
grid on

Aw = linspace(0,80,200);
t_fixed = 0.15;

figure
plot(Aw, C(t_fixed,Aw), 'LineWidth',2)
xlabel('Window area A_w (m^2)')
ylabel('Cost C (€)')
title('Cost vs Window Area (t = 0.15 m)')
grid on
