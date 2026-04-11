function [x_best,cost_best,T_current] = simulated_annealing(f, x0, T0, Tmin, alpha, maxIter, epsilon, g, h,lambda0,variablePenalty)

    % Initialization
    x_current = x0;
    T_current = T0;

    

    cost_current = f(x_current) + (lambda0/T_current)*penalty(x_current, g, h);

    x_best = x_current;
    cost_best = cost_current;

    % Main loop
    while T_current > Tmin
        for k = 1 : maxIter

            % Generate random neighbor
            x_new = x_current;
            sz = size(x_new);
            j = randi([1, sz(1)]);

            delta_x = (2 * rand() - 1) * epsilon;
            x_new(j) = x_new(j) + delta_x;

            % New cost WITH penalty
            if variablePenalty
                lambda = lambda0 / T_current;
            else
                lambda = lambda0;
            end
            cost_new = f(x_new) + lambda * penalty(x_new, g, h);

            delta_cost = cost_new - cost_current;

            if delta_cost <= 0
                x_current = x_new;
                cost_current = cost_new;
            else
                p = exp(-delta_cost / T_current);
                if rand() < p
                    x_current = x_new;
                    cost_current = cost_new;
                end
            end

            % Best solution
            if cost_current < cost_best
                x_best = x_current;
                cost_best = cost_current;
            end
        end

        T_current = alpha * T_current;
    end


    % Penalty function
    function P = penalty(x, g, h)
        P = 0;

        % Inequality constraints
        for i = 1:length(g)
            gi = g{i}(x);
            P = P + max(0, gi)^2;
        end

        % Equality constraints
        for j = 1:length(h)
            hj = h{j}(x);
            P = P + hj^2;
        end
    end

end