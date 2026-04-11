while true

    choice = menu("Select an option", ...
        'Initial Problem Investigation', ...
        'Validation on Himmelblau Function', ...
        'Optimization of Simplified Problem', ...
        'Optimization of Full Problem', ...
        'Exit');

    switch choice
        case 1
            initial_problem_investigation

        case 2
            initial_validation

        case 3
            simplified_problem

        case 4
            full_problem

        case 5
            disp("Exiting program...");
            break

        otherwise
            disp("Invalid Choice")
    end

end