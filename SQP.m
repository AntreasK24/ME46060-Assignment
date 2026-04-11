function [x_best,cost_best] = SQP(f,gradf,x0,confun,jaccon,debug)
    
    %Initialize
    x = x0;
    n = length(x);
    B = eye(n);

    rho = 1000;

    tol_kkt = 1e-6;

    options_qp = optimoptions('quadprog','Display','off');


    cost = f(x);

    [Jineq, Jeq] = jaccon(x);
    [cineq,ceq] = confun(x);

    lambda = zeros(length(ceq),1);
    mu = zeros(length(cineq),1);



    max_iter = 100;
    tol = 1e-6;


    phi = @(x) merit(x);

    function val = merit(x)
        [cineq, ceq] = confun(x);
        val = f(x) + rho * sum(max(0, cineq)) + rho * sum(abs(ceq));
    end

    function cin = get_cineq(x)
        [cin, ~] = confun(x);
    end

    %QP Step
    iter = 1;
    while iter < max_iter
        alpha = 1;

        % 1. Evaluate at current x
        g = gradL(x, lambda, mu);
        [cineq, ceq] = confun(x);
        [Jineq, Jeq] = jaccon(x);

        % 2. Build QP
        H = B;
        % Build QP constraints safely

        if isempty(cineq)
            A = [];
            b = [];
        else
            A = Jineq;
            b = -cineq + 1e-6; % small relaxation
        end

        if isempty(ceq)
            Aeq = [];
            beq = [];
        else
            Aeq = Jeq;
            beq = -ceq;
        end

        % 3. Solve QP
        [d,fval,exitflag,output,lambda_QP] = quadprog(H, g, A, b, Aeq, beq, [], [], [], options_qp);
        if debug
            fprintf('QP exitflag: %d\n', exitflag);
        end

        % --- Multiplier handling (ROBUST FIX) ---
        if exitflag ~= 1 || isempty(lambda_QP)
            % QP failed → fallback
            warning('QP failed → using zero multipliers');

            mu_new = zeros(size(cineq));
            lambda_new = zeros(size(ceq));

        else
            if isstruct(lambda_QP)
                mu_new = lambda_QP.ineqlin;
                lambda_new = lambda_QP.eqlin;
            else
                % fallback for old MATLAB (rare case)
                mu_new = zeros(size(cineq));
                lambda_new = zeros(size(ceq));
            end
        end




        if norm(d) < 1e-8
            disp('Step too small → stopping (likely stagnation)');
            break;
        end

        alpha = 1;
        phi_current = phi(x);

        while true
            x_trial = x + alpha*d;

            if phi(x_trial) < phi_current
                break;
            end

            alpha = alpha/2;

            if alpha < 1e-6
                break;
            end
        end

        x_new = x + alpha*d;


        s = x_new - x;

        gL_old = gradL(x, lambda, mu);
        gL_new = gradL(x_new, lambda_new, mu_new);

        y = gL_new - gL_old;

        if norm(s) > 1e-8 && s'*y > 1e-10 && s'*B*s > 1e-12
            B = B ...
                - (B*s*s'*B)/(s'*B*s) ...
                + (y*y')/(s'*y);

            B = 0.5*(B + B');
        end

        x = x_new;

        lambda = lambda_new;
        mu = mu_new;

        cost = f(x);


        iter = iter + 1;

        %KKT Criteria
        gL = gradL(x,lambda,mu);
        r_stationarity = norm(gL);


        if norm(d) < 1e-8 && r_stationarity > 1e-3
            disp('Stuck → resetting B and continuing');
            B = eye(n);
            iter = iter + 1;
            continue;
        end

        [cineq, ceq] = confun(x);

        r_eq = norm(ceq);
        r_ineq = norm(max(0, cineq));

        r_comp = norm(mu .* cineq);

        r_dual = norm(min(0, mu));

        fprintf('Iter %d: stat=%.2e, eq=%.2e, ineq=%.2e, comp=%.2e\n', ...
        iter, r_stationarity, r_eq, r_ineq, r_comp);

        %Stop if criteria is satisfied
        if r_stationarity < tol_kkt && ...
            r_eq < tol_kkt && ...
            r_ineq < tol_kkt && ...
            r_comp < tol_kkt && ...
            r_dual < tol_kkt
            break;
        end



    end
    x_best = x;
    cost_best = cost;

    function gL = gradL(x, lambda, mu)
        gL = gradf(x);

        [Jineq, Jeq] = jaccon(x);

        % Add equality part ONLY if it exists
        if ~isempty(Jeq) && ~isempty(lambda)
            gL = gL + Jeq' * lambda;
        end

        % Add inequality part ONLY if it exists
        if ~isempty(Jineq) && ~isempty(mu)
            gL = gL + Jineq' * mu;
        end
    end


end

