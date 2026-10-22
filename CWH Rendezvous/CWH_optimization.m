%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA Optimization Intro
% Convex Optimization in ASA Example
% Author: Travis Hastreiter 
% Created On: 17 October, 2026
% Description: CWH Convex Trajectory Optimization with ZOH control
% Define the target spacecraft, the ISS, as the Chief
% Define the docking spacecraft that we control as the Deputy
% Assumptions:
% - Deputy is in a circular orbit
% - No external perturbations (except thrust for deputy)
% - No mass change from thrust so constant mass
% Most Recent Change: 21 October, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% Physical Parameters
mu = 398600.4415; % [km2 / s3] Earth gravitational parameter

% Boundary Conditions
r_0 = [-1.5; -0.5; 0.2]; % [km]
v_0 = [3e-3; 0; 0]; % [km / s]

r_f = [0; 0; 0]; % [km]
v_f = [0; 0; 0]; % [km / s]

x_0 = [r_0; v_0];
x_f = [r_f; v_f];

% Timespan
t_0 = 0; % [s]
t_f = 3600; % [s]

t_span = [t_0, t_f];

% Constraint parameters
r_max = 3.5; % [km]
u_max = 5e-6; % [km / s^2] control is acceleration from thrust

% Problem Parameters
N = 50; % Number of time discretization nodes
Nu = N - 1; % Number of control nodes

nx = 6; % [r_r; r_t; r_n; v_r; v_t; v_n] state vector dimension
nu = 3; % [u_r; u_t; u_n] control vector dimension
np = 0; % [] parameter vector dimension - (no parameters, such as final time, for this problem)

t_k = linspace(t_span(1), t_span(2), N); % Discretize time into intervals
delta_t = t_k(2) - t_k(1); % Get the constant time interval size

% Algorithm Parameters
default_tolerance = 1e-12; % Integration tolerance
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% Chief Orbit Parameters
r_c = 6728; % [km] ISS Keplerian circular orbit

n = sqrt(mu / r_c ^ 3); % Calculate mean motion of orbit (equal to angular velocity here since it is a circular orbit)

%% Linearization
f = @(t, x, u, p) CWH_EoM(t, x, u, mu, r_c); % Equations of motion (can write as f = Ax + Bu here since it is linear)
A = @(t, x, u, p) [zeros(3), eye(3); [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2], [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0]]; % Manual partial f partial x jacobian
B = @(t, x, u, p) [zeros(nu); eye(nu)]; % Manual partial f partial u jacobian
E = @(t, x, u, p) zeros([nx, np]); % Manual partial f partial p jacobian
c = @(t, x, u, p) f(t, x, u, p) - (A(t, x, u, p) * x + B(t, x, u, p) * u + E(t, x, u, p) * p); % Always going to be zero because system is linear
    
% Create reference trajectory (can be anything because dynamics are linear and constraints are convex
x_ref = ones([nx, N]);
u_ref = zeros([nu, Nu]);
p_ref = zeros([np, 1]); % Parameters are constants and do not depend on time

%% Discretization
[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_span, x_ref, u_ref, p_ref, tolerances);

%% Solve with CVX
cvx_begin
    variable X_sol(nx, N)
    variable U_sol(nu, N - 1)
    minimize( sum(norms(U_sol, 2, 1) * delta_t) )
    subject to
        % Dynamics constraints
        for k = 1 : (N - 1)
            X_sol(:, k + 1) == A_k(:, :, k) * X_sol(:, k) + B_k(:, :, k) * U_sol(:, k) + c_k(:, :, k);
        end

        % Path constraints 
        norms(X_sol(1:3, :), 2, 1) <= r_max; % Max distance constraint
        norms(U_sol, 2, 1) <= u_max; % Max control magnitude constraint

        % State Boundary Conditions
        X_sol(:, 1) == x_0; % Initial condition
        X_sol(:, end) == x_f; % Terminal condition
cvx_end

%% Plot Solution
figure
plot3(X_sol(1, :), X_sol(2, :), X_sol(3, :)); hold on
quiver3(X_sol(1, 1:Nu), X_sol(2, 1:Nu), X_sol(3, 1:Nu), U_sol(1, :), U_sol(2, :), U_sol(3, :)); hold off
title("Convex Minimum Fuel Control of Relative Orbit Transfer Near ISS")
xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
legend("Trajectory", "Thrust", location = "eastoutside")
grid on

figure
tiledlayout(1, 3, "TileSpacing","compact")

nexttile
plot(t_k, X_sol(1, :), DisplayName="r"); hold on
plot(t_k, X_sol(2, :), DisplayName="\theta"); hold on
plot(t_k, X_sol(3, :), DisplayName="n"); hold on
r_mag = vecnorm(X_sol(1:3, :), 2, 1);
plot(t_k, r_mag, DisplayName="Distance"); hold on
yline(r_max, DisplayName="Max Distance")
title("RTN Frame Position History")
xlabel("Time [s]")
ylabel("Position [km]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

nexttile
plot(t_k, X_sol(4, :), DisplayName="rdot"); hold on
plot(t_k, X_sol(5, :), DisplayName="\theta dot"); hold on
plot(t_k, X_sol(6, :), DisplayName="ndot"); hold off
title("RTN Frame Velocity History")
xlabel("Time [s]")
ylabel("Velocity [km / s]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

nexttile

u_mag = vecnorm(U_sol, 2, 1);

stairs(t_k(1:Nu), U_sol(1, :) * 1e6, DisplayName="u_r"); hold on
stairs(t_k(1:Nu), U_sol(2, :) * 1e6, DisplayName="u_\theta"); hold on
stairs(t_k(1:Nu), U_sol(3, :) * 1e6, DisplayName="u_n"); hold on
stairs(t_k(1:Nu), u_mag * 1e6, DisplayName="||u||"); hold off
title("RTN Frame Control History")
xlabel("Time [s]")
ylabel("Control [mm / s2]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

sgtitle("State and Control Time Histories for Relative Orbit Transfer Around ISS")

%% Helper Functions
function [xdot] = CWH_EoM(t, x, u, mu, r_c)
    r = x(1:3);
    v = x(4:6);

    n = sqrt(mu / r_c ^ 3);

    rdot = v;
    vdot = [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2]  * r ...
        +  [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0] * v ...
        + u;

    xdot = [rdot; vdot];
end