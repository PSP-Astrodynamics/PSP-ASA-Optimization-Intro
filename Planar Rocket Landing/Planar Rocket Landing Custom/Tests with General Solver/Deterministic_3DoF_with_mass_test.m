%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA Optimization Intro
% SCP Planar Rocket Landing Optimization
% Author: Travis Hastreiter
% Created On: 13 November, 2026
% Description: 3DoF planar rocket landing with changing mass
% Most Recent Change: 14 November, 2026
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
 
% The vehicle object is only created because the initial guess function
% used wants it
vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I, alpha = alpha);

% Problem Parameters
tf = 6; % [s]
N = 30; % []
r_0 = [-10; 40]; % [km]
theta_0 = deg2rad(100); % [rad]
v_0 = [2; -5]; % [km / s]
w_0 = 0.01; % [rad / s]
glideslope_angle_max = deg2rad(45); % [rad]

nx = 7;
nu = 2;
np = 0;

x_0 = [r_0; v_0; theta_0; w_0; m_0];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

initial_guess = "straight line"; % "CasADi" or "straight line" or "screw"
parser = "CVX";

% PTR algorithm parameters
ptr_ops.iter_max = 15;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 1e-4;
ptr_ops.w_vc = 1e4;
ptr_ops.w_tr = ones(1, Nu) * 5e-1;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 2e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 0;
ptr_ops.alpha_p = 0;

scale = false;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsPlanar3DoF(t, x, u, L, I, Isp);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = {1:N, @(t, x, u, p) norm(x(1:2)) - x(2) / cos(glideslope_angle_max)};
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = {1:N, @(t, x, u, p) norm(u) - T_max};
max_gimbal_constraint = {1:N, @(t, x, u, p) norm(u) - u(1) / cos(gimbal_max)};
control_convex_constraints = {max_gimbal_constraint,max_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
glideslope_STC_angle = deg2rad(5);
glideslope_STC_trigger_height = 0.5; % [km]
glideslope_STC_constraint = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) (norm(x(1:2, k)) - x(2, k) / cos(glideslope_STC_angle)) * (x_ref(2, k) <= glideslope_STC_trigger_height)};
state_nonconvex_constraints = {};

nonconvex_constraints = [state_nonconvex_constraints];

% Terminal boundary conditions
terminal_bc = @(x, p, x_ref, p_ref) [x(1:6, :) - x_f; 0];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(norms(u, 2, 1)) * delta_t / T_max;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * delta_t;
end

%% Create Guess
sl_guess = guess_3DoF(x_0(1:6), x_f + [0; 0; 0; 0; 0; 0], N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [0, cumsum(sl_guess.u(3, 1:(end - 1)) * delta_t), sum(sl_guess.u(3, 1:(end - 1))) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [0, cumsum(sl_guess.u(3, 1:(end - 1)) * delta_t)]];
end
sl_guess.u = sl_guess.u(1:nu, :);

if initial_guess == "straight line"
    guess = sl_guess;
elseif initial_guess == "CasADi"
    guess = CasADi_sol;
end
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;
% 
% figure
% plot_3DoFc_trajectory(t_k, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_3DoFc_time_histories(t_k, sl_guess.x, sl_guess.u)
% 
% figure
% plot_3DoFc_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)
% 
% figure
% plot_3DoFc_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc, nonconvex_constraints = nonconvex_constraints, integration_tolerance = 1e-8, discretization_method = "error", N_sub = 1);

%% Test Discretization on Initial Guess

[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_3DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.u, guess.p);
% 
% figure
% comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")
% 
% figure
% comparison_plot_3DoFc_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%%
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, nonconvex_constraints = nonconvex_constraints, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1);
prob_3DoF.vehicle = vehicle;

[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

%% Solve Problem with PTR
ptr_sol_vc = ptr(prob_3DoF, ptr_ops, parser);

%%
ptr_sol = ptr_sol_vc;

if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

%%
figure
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_3DoF.objective(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p), [ptr_sol.info.J]]); hold on
%yline(CasADi_sol.objective); hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
yscale("log")
title("Stopping Criteria vs Iteration")
grid on

nexttile
plot(0:ptr_sol.converged_i, squeeze(sum(vecnorm(ptr_sol.Delta(:, :, 1:(ptr_sol.converged_i + 1)), 2, 1))))
yscale("log")
title("Defect Norm vs Iteration")
grid on

% %%
% Js = [ptr_sol.info.J];
% (Js(end) - CasADi_sol.objective) / CasADi_sol.objective * 100

%%
figure
plot(0:ptr_sol.converged_i, [prob_3DoF.objective(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p), [ptr_sol.info.J]]); hold on
%yline(CasADi_sol.objective); hold off
xlabel("PTR Iteration")
ylabel("\Delta V [km / s]")
legend("PTR Iterations", "CasADi Solution", location = "southeast")
grid on

%%
figure
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
xlabel("PTR Iteration")
ylabel("Defect")
yscale("log")
xlim([0, ptr_sol.converged_i])
grid on

%%
i = ptr_sol.converged_i;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p);

plot_3DoFc_trajectory(t_cont_sol, x_cont_sol, u_cont_sol, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

%%

figure
plot_3DoFc_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)
% %%
% figure
% comparison_plot_3DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i), CasADi_sol.x}, ["Guess", "Continuous Propagation", "Solution Output", "CasADi"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "3DoF Solution Comparison")
% %%
figure
comparison_plot_3DoF_trajectory({guess.x, x_cont_sol}, ["Guess", "PTR"], glideslope_angle_max, linestyle = [":", "-"], title = "", STC_glideslope_angle_height = [glideslope_STC_angle, glideslope_STC_trigger_height])

%%
figure
comparison_plot_3DoFc_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

% %%
% comparison_plot_3DoFc_time_histories({t_k, t_cont_sol, t_k}, {CasADi_sol.x, x_cont_sol, ptr_sol.x(:, :, i)}, {CasADi_sol.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["CasADi", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

%%
t_iters = {};
x_iters = {};
u_iters = {};
linestyle = string(1:ptr_sol.converged_i);
for i = 1:ptr_sol.converged_i
    t_iters{i} = t_k;
    x_iters{i} = ptr_sol.x(:, :, i);
    u_iters{i} = ptr_sol.u(:, :, i);
    linestyle(i) = "-";
end
figure
comparison_plot_3DoF_trajectory(x_iters, "iter " + string(1:ptr_sol.converged_i), glideslope_angle_max, linestyle = linestyle, title = "Solution vs Iteration")

figure 
comparison_plot_3DoFc_time_histories(t_iters, x_iters, u_iters, "iter " + string(1:ptr_sol.converged_i), linestyle = linestyle, title = "Solution vs Iteration")
