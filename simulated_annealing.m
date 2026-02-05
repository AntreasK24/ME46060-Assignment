function [x_best,cost_best,T_current] = simulated_annealing(f, x0, T0, Tmin, alpha, maxIter, epsilon)

    % Initialization
    x_current = x0;
    T_current = T0;

    cost_current = f(x_current);

    x_best = x_current;
    cost_best = cost_current;

    % Main loop
    while T_current > Tmin
        for k = 1 : maxIter
            %Generate random neighbor
            x_new = x_current;
            h = size(x_new);
            j = randi([1, h(1)]);

            delta_x = (2 * rand() - 1) * epsilon;


            x_new(j) = x_new(j) + delta_x;

            cost_new = f(x_new);

            delta_cost = cost_new - cost_current;

            if delta_cost <= 0
                x_current = x_new;
                cost_current = cost_new;
            else
                p = exp( -delta_cost / T_current );
                if rand([0,1]) < p
                    x_current = x_new;
                    cost_current = cost_new;
                end
            end



            %Best solution

            if cost_current < cost_best
                x_best = x_current;
                cost_best = cost_current;
            end
        end

        T_current = alpha * T_current;


    end

end
