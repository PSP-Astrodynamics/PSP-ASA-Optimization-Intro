%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA Optimization Intro
% Planar 3DoF Rocket Landing (without general SCP solver)
% Author: Travis Hastreiter 
% Created On: 12 November, 2026
% Description: Planar 3DoF rocket landing using specialized SCP solver
% Most Recent Change: 12 November, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Vehicle Parameters
Isp = 200; % [s]
alpha = 1 / (9.81 * Isp); % [s / km]
m_dry = 100; % [kg]
m_wet = 200; % [kg]
m_0 = m_wet;
T_max = 2500; % [kg km / s2]
T_min = 0.55 * T_max; % [kg km / s2]
R = 0.5; % [m]
I = m_0 * R ^ 2 / 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 0.2; % [km] Distance from CoM to nozzle
gimbal_max = deg2rad(10); % [rad]
 
% Initial condition
r_0 = [-10; 40]; % [km]
theta_0 = deg2rad(100); % [rad]
v_0 = [2; -5]; % [km / s]
w_0 = 0.01; % [rad / s]

x_0 = [r_0; v_0; theta_0; w_0; m_0];

% Terminal condition
r_f = [0; 0]; % [km]
theta_f = deg2rad(90); % [rad]
v_f = [0; 0]; % [km / s]
w_f = 0; % [rad / s]

x_f = [r_f; v_f; theta_f; w_f];

% Problem Parameters
tf = 60; % [s]
N = 30; % []
glideslope_max = deg2rad(45); % [rad]

nx = 7;
nu = 2;
np = 0;

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
t_cont = linspace(tspan(1), tspan(2), 200);
delta_t = t_k(2) - t_k(1);
integration_tolerance = 1e-12;
tolerances = odeset(RelTol = integration_tolerance, AbsTol = integration_tolerance);

Nu = N - 1;

% PTR algorithm parameters
ptr_ops.iter_max = 25;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 4e-5;
ptr_ops.w_vc = 1e4;
ptr_ops.w_tr = 1e-1;
ptr_ops.delta_tol = 5e-3;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1e-6;

%% Specify Objective
objective = @(x, u) sum(norms(u, 2, 1) * delta_t) / T_max;

%% Create Initial Guess
% Linear interpolation
r_guess = (x_f(1:2) - x_0(1:2)) / tf * t_k + x_0(1:2);
v_guess = (x_f(1:2) - x_0(1:2)) / tf .* ones(2, N);
theta_guess = (x_f(5) - x_0(5)) / tf * t_k + x_0(5);
w_guess = (x_f(5) - x_0(5)) / tf .* ones(1, N);
m_guess = linspace(m_dry, m_wet, N);

% Assemble initial guess
guess.x = [r_guess; v_guess; theta_guess; w_guess; m_guess];
guess.u = zeros([nu, Nu]) + [T_max/2; 1e-6];
guess.p = zeros([np, 1]); % Parameters are constants and do not depend on time

% Plot initial guess
figure
plot_3DoF_trajectory(t_k, guess.x, guess.u, glideslope_max, gimbal_max, 0, T_max)
figure
plot_3DoF_time_histories(t_k, guess.x, guess.u)

%% Load dynamics
f = @(t, x, u, p) SymDynamicsPlanar3DoF(t, x, u, L, I, Isp);
A = @(t, x, u, p) SymAJacobianPlanar3DoF(t, x, u, L, I, Isp);
B = @(t, x, u, p) SymBJacobianPlanar3DoF(t, x, u, L, I, Isp);
E = @(t, x, u, p) zeros([nx, np]); % Manual partial f partial p jacobian
c = @(t, x, u, p) f(t, x, u, p) - (A(t, x, u, p) * x + B(t, x, u, p) * u + E(t, x, u, p) * p); % Always going to be zero because system is linear

%% Test Out Discretization
[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, tspan, guess.x, guess.u, guess.p, tolerances);

% Propagate initial guess with discretized dynamics
x_disc = disc_prop(x_0, guess.u, guess.p, A_k, B_k, E_k, c_k);

% Propagate initial guess with continuous dynamics
[x_cont_guess, u_cont_guess] = cont_prop(f, t_k, t_cont, x_0, guess.u, guess.p, tolerances);

% Plot initial guess vs propagated initial guess
figure
comparison_plot_3DoF_trajectory({guess.x, x_cont_guess, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Sequential Convex Programaming using Penalized Trust Region (PTR) method
% Initialize outputs
delta_x = zeros(1, ptr_ops.iter_max);
total_defects = zeros(1, ptr_ops.iter_max + 1);
objectives = zeros(1, ptr_ops.iter_max);

% Set initial reference as the initial guess
x_ref = guess.x;
u_ref = guess.u;
p_ref = guess.p;
[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, tspan, x_ref, u_ref, p_ref, tolerances);
total_defects(1) = sum(norms(Delta, 2, 1));

% PTR loop
for i = 1 : ptr_ops.iter_max
    % Solve convex subproblem
    [X_sol, U_sol, V_sol] = solve_ptr_subproblem(objective, x_ref, u_ref, A_k, B_k, c_k, x_0, x_f, glideslope_max, gimbal_max, T_max, N, delta_t, ptr_ops);

    delta_x(i) = sum(sum_square(X_sol - x_ref)); % Change from previous solution
    objectives(i) = objective(X_sol, U_sol);

    % Linearize and discretize around new solution
    x_ref = X_sol;
    u_ref = U_sol;
    [A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, tspan, x_ref, u_ref, p_ref, tolerances);
    
    % Sum total defects (difference between predicted and integrated final state) over each time interval
    total_defects(i + 1) = sum(norms(Delta, 2, 1));

    % Print iteration results
    fprintf("PTR Iter %.2g - J: %.5g, Delta: %.3g, delta_x: %.3g\n", i, objectives(i), total_defects(i + 1), delta_x(i))
    
    % End SCP loop if tolerance is met
    if total_defects(i + 1) < ptr_ops.Delta_min && delta_x(i) < ptr_ops.delta_tol && i >= ptr_ops.iter_min
        break;
    end
end

%% Plot PTR Convergence
figure
tiledlayout(1, 3)

nexttile
plot(1 : i, objectives(1 : i)); hold on
legend("PTR Iterations")
title("Objective vs Iteration")
grid on

nexttile
plot(delta_x(1 : i))
yscale("log")
title("Change in State Trajectory vs Iteration")
grid on

nexttile
plot(0 : i, total_defects(1 : (i + 1)))
yscale("log")
title("Defect Norm vs Iteration")
grid on

%% Plot Solution
[x_cont, u_cont] = cont_prop(f, t_k, t_cont, x_0, U_sol, guess.p, tolerances);

figure
plot_3DoF_trajectory(t_cont, x_cont, u_cont, glideslope_max, gimbal_max, 0, T_max)
figure
plot_3DoF_time_histories(t_cont, x_cont, u_cont)

%% Helper Functions
function [x_cont, u_cont] = cont_prop(f, t_k, t_cont, x_0, u_k, p, tolerances)
    Nu = size(u_k, 2);
    u_func = @(t) interp1(t_k(1:Nu), u_k', t, "previous", "extrap")';

    [~, x_cont] = ode45(@(t, x) f(t, x, u_func(t), p), t_cont, x_0, tolerances);
    x_cont = x_cont';

    u_cont = u_func(t_cont(1:(numel(t_cont) - 1)));
end

function [x_disc] = disc_prop(x_0, u, p, A_k, B_k, E_k, c_k)
    nx = numel(x_0);
    N = size(A_k, 3) + 1;

    x_disc = zeros(nx, N);
    x_disc(:, 1) = x_0;
    
    for k = 1:(N - 1)
        x_disc(:, k + 1) = A_k(:, :, k) * x_disc(:, k) ...
                     + B_k(:, :, k) * u(:, k) ...
                     + zero_if_empty(E_k(:, :, k) * p) ...
                     + c_k(:, :, k);
    end
end

function [X_sol, U_sol, V_sol] = solve_ptr_subproblem(objective, x_ref, u_ref, A_k, B_k, c_k, x_0, x_f, glideslope_max, gimbal_max, T_max, N, delta_t, ptr_ops)
    nx = size(x_ref, 1);
    nu = size(u_ref, 1);
    Nu = N - 1;

    cvx_begin quiet
    variable X_sol(nx, N)
    variable U_sol(nu, Nu)
    variable V_sol(nx, N - 1)
    minimize( objective(X_sol, U_sol) ... % Objective (actual function we want to minimize)
            + sum(norms(V_sol, 1, 1)) * ptr_ops.w_vc ... % Virtual control cost
            + sum(sum_square(X_sol - x_ref)) * ptr_ops.w_tr * ptr_ops.alpha_x ... % State trust region cost
            + sum(sum_square(U_sol - u_ref)) * ptr_ops.w_tr * ptr_ops.alpha_u) % Control trust region cost
    subject to
        % Dynamics constraints
        for k = 1 : (N - 1)
            X_sol(:, k + 1) == A_k(:, :, k) * X_sol(:, k) + B_k(:, :, k) * U_sol(:, k) + c_k(:, :, k) + V_sol(:, k);
        end

        % Convex State Path Constraints 
        norms(X_sol(1:2, :), 2, 1) <= X_sol(2, :) / cos(glideslope_max); % Glideslope constraint
        
        % Convex Control Path Constraints
        norms(U_sol, 2, 1) <= T_max; % Max thrust constraint
        norms(U_sol(1:2, :), 2, 1) <= U_sol(1, :) / cos(gimbal_max); % Gimbal constraint

        % State Boundary Conditions
        X_sol(:, 1) == x_0; % Initial condition
        X_sol(1:(nx - 1), end) == x_f; % Terminal condition (final mass unconstrained)
    cvx_end
end