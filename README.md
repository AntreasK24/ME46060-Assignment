# ME46060-Assignment

This repository holds all the scripts necessary to solve an optimization problem for a house's heating system. For details on the problem formulation and optimization choices please refer to the corresponding report. 

## How to run

In order to execute the problem run the `main.m` script this will present you with a menu with the following options

### Iniitial Problem Investigation
Checks the problems monotonicity, convexity, boundedness and sensitivity numerically by running `initial_problem_investigation.m`.

### Validation on Himmelblau Function
Runs our proposed optimization approach on Himmelblau's function as a way of validating our approach by running `initial_validation.m`.

### Optimization of simplified problem
Runs our proposed optimization approach on a simplified version of the problem with only two decision variables. This is to further validate our approach and be able to visualize the results. This option runs `simplified_problem.m`

### Optimization of full problem
Runs our proposed optimization approach on the full version of the problem with 5 decision variables and checks if the KKT conditions are satisfied. This option runs `full_problem.m`

### Exit
Exits the script.