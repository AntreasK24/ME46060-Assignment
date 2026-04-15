function [x_best,cost_best] = SQP(f,gradf,x0,confun,jaccon,debug,rho)
    
    %Initialize
    x = x0;
    n = length(x);
    B = eye(n);

    

    % options_qp = optimoptions('quadprog','Display','off');
    options_qp = optimset('Display','off');


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

        % Evaluate at current x
        g = gradL(x, lambda, mu);
        [cineq, ceq] = confun(x);
        [Jineq, Jeq] = jaccon(x);

        % Build QP
        H = B;
        % Build QP constraints safely

        if isempty(cineq)
            A = [];
            b = [];
        else
            A = Jineq;
            b = -cineq;
        end

        if isempty(ceq)
            Aeq = [];
            beq = [];
        else
            Aeq = Jeq;
            beq = -ceq;
        end

        % Solve QP
        [d,fval,exitflag,output,lambda_QP] = quadprog(H, g, A, b, Aeq, beq, [], [], [], options_qp);
        if debug
            fprintf('QP exitflag: %d\n', exitflag);
        end

        % Multiplier handling
        if exitflag ~= 1 || isempty(lambda_QP)
            warning('QP failed → using zero multipliers');

            mu_new = zeros(size(cineq));
            lambda_new = zeros(size(ceq));

        else
            if isstruct(lambda_QP)
                mu_new = lambda_QP.ineqlin;
                lambda_new = lambda_QP.eqlin;
            else
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

        viol = @(x) norm(max(0, get_cineq(x)),1);

        while true
            x_trial = x + alpha*d;

            if phi(x_trial) < phi_current && viol(x_trial) <= viol(x)
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




        if norm(d) < 1e-8 && r_stationarity > 1e-3
            disp('Stuck → resetting B and continuing');
            B = eye(n);
            iter = iter + 1;
            continue;
        end

        [cineq, ceq] = confun(x);


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

