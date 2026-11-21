%% Initialize
% Mass
m_dry = 100; % [kg]
m_wet = 200; % [kg]

% Thrust
T_max = 2500; % [N]
T_min = 0.55*T_max; % [N]
isp = 200; % [sec]

% Rotation 
R = 0.5; %[m]
I = (m_wet^2) / R; %[kg*m^2]
L = 0.2; %[m]

%% Inital Conditions
r_0 = [-10; 40]; % [m]
v_0 = [2; -5]; % [m/s]
theta_0 = deg2rad(100); % [rad]
w_0 = 0.01; % [rad / s]
x_0 = [r_0; v_0; theta_0; w_0; m_wet]; 

%% Terminal Conditions
r_f = [0; 0]; % [m]
v_f = [0; 0]; % [m/s]
theta_f = deg2rad(90); % [rad]
w_f = 0; % [rad / s]
x_f = [r_f; v_f; theta_f; w_f; m_dry];

%% Problem Parameters
t_f = 7; % [s]
t_0 = 0; % [s]
t_k = linspace(t_0, t_f, N); 
delta_t = t_k(2) - t_k(1); % [s]

% Dimensions
nx = 7;
nu = 2;
np = 0;

% Angles
glide_slope = deg2rad(45);
gimbal_angle = deg2rad(10);

% Create Guesses/References
N = 30; 
Nu = N - 1;
r_guess = (r_f - r_0) / t_f * t_k + r_0; 
v_guess = ((r_f - r_0) / t_f) .*ones(2, N);

theta_guess = ((theta_f - theta_0) / t_f) * t_k + theta_0;
w_guess = ((theta_f - theta_0) / t_f) .*ones(1, N);
m_guess = ((m_dry - m_wet) / t_f) * t_k + m_wet;

x_guess = [r_guess; v_guess; theta_guess; w_guess; m_guess];
u_guess = [9.81*298; 1] .*ones(2, N-1);
p_guess = zeros(0, 1);

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsPlanar3DoF(t, x, u, L, I, Isp);
A = @(t, x, u, p) SymAJacobianPlanar3DoF(t, x, u, L, I, Isp);
B = @(t, x, u, p) SymBJacobianPlanar3DoF(t, x, u, L, I, Isp);
E = @(t, x, u, p) zeros([nx, np]); % Manual partial f partial p jacobian
c = @(t, x, u, p) f(t, x, u, p) - (A(t, x, u, p) * x + B(t, x, u, p) * u + E(t, x, u, p) * p); % Always going to be zero because system is linear

% f = @(t, x, u, p) dynamics(t, x, u, L, I, isp);
% A = @(t, x, u, p) Jacobian_x(t, x, u, L, I, isp);
% B = @(t, x, u, p) Jacobian_u(t, x, u, L, I, isp);
% E = @(t, x, u, p) zeros([nx, np]);
% c = @(t, x, u, p) f(t, x, u, p) - (A(t, x, u, p) * x + B(t, x, u, p) * u + E(t, x, u, p) * p);

%% Discretization
integration_tolerance = 1e-12;

tolerances = odeset(RelTol=integration_tolerance, AbsTol=integration_tolerance);
[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_k, x_guess, u_guess, p_guess, tolerances);

%% Integration
u_func = @(t) interp1(t_k(1:Nu), u_guess', t, "previous", "extrap")';
[t_cont, x_cont] = ode45(@(t, x) f(t, x, u_func(t), p_guess), t_k, x_0, tolerances);
x_cont = x_cont';

figure
plot(x_cont(1, :), x_cont(2, :));
grid on;
axis equal;
title('Integrated Intial Guess');
xlabel('X [m]');
ylabel('Y [m]');

%% Solve with CVX
[X_sol, U_sol] = solving(1e-1, 1e-7, 1e4, 0.001, 0.001, 25, Delta, nx, nu, N, delta_t, x_guess, u_guess, glide_slope, T_max, gimbal_angle, m_dry, x_0, x_f, f, A, B, E, c, t_k, p_guess, tolerances);


%% Graphing
y_max = max(X_sol(2, :));
x_lim = [min(X_sol(1, :)) - y_max / 10, max(X_sol(1, :)) + y_max / 10];
y_lim = [min(X_sol(2, :)) - y_max / 10, y_max * 1.1];

options.step = 4;
options.title = "State Trajectory";

figure
plot(X_sol(1, :), X_sol(2, :)); hold on
line([x_lim(1), 0, x_lim(2)], abs([x_lim(1), 0, x_lim(2)]) / tan(glide_slope), 'Color', 'w', 'LineStyle', '--'); hold on
quiver(X_sol(1, 1:Nu), X_sol(2, 1:Nu), -U_sol(1, 1:end), -U_sol(2, 1:end), ShowArrowHead = "off", color = "r", AutoScaleFactor=0.4)
title(options.title)
legend("", "Glideslope", "Thrust", Location="southoutside", Orientation="horizontal")
xlabel("X [km]")
ylabel("Y [km]")
axis equal
grid on
xlim(x_lim)
ylim(y_lim)

%% Graphing Position, Velocity, Orientation, Angular Velocity, Mass vs Time
figure
tiledlayout(4, 1);
sgtitle('Position, Velocity, Orientation, Angular Velocity vs Time');

% Position vs. Time
nexttile;
plot(t_k, X_sol(1, :), 'LineWidth', 1.5); hold on;
plot(t_k, X_sol(2, :), 'LineWidth', 1.5);
grid on;
title('Position vs. Time');
xlabel('Time [s]');
ylabel('Position [m]');
legend('X Position', 'Y Position', 'Location', 'best');

% Velocity vs. Time
nexttile;
plot(t_k, X_sol(3, :), 'LineWidth', 1.5); hold on;
plot(t_k, X_sol(4, :), 'LineWidth', 1.5);
grid on;
title('Velocity vs. Time');
xlabel('Time [s]');
ylabel('Velocity [m]');
legend('X Velocity', 'Y Velocity', 'Location', 'best');

% Orientation vs. Time
nexttile;
plot(t_k, X_sol(5, :), 'LineWidth', 1.5);
grid on;
title('Orientation vs. Time');
xlabel('Time [s]');
ylabel('Orientation [rad]');

% Angular Velocity vs. Time
nexttile;
plot(t_k, X_sol(5, :), 'LineWidth', 1.5);
grid on;
title('Angular Velocity vs. Time');
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');

% Mass vs. Time
figure
plot(t_k, X_sol(7, :), 'LineWidth', 1.5);
grid on;
title('Mass vs. Time');
xlabel('Time [s]');
ylabel('Mass [kg]');
%%
% Thrust, Gimbal vs. Time
figure
tiledlayout(2, 1);
sgtitle('Thrust, Gimbal vs. Time')

nexttile;
plot(t_k(1:Nu), norms(U_sol, 2, 1), 'LineWidth', 1.5); hold on;
grid on;
title('Thrust vs. Time');
xlabel('Time [s]');
ylabel('Thrust [N]');


nexttile;
plot(t_k(1:Nu), atan2(U_sol(2, :), U_sol(1, :)), 'LineWidth', 1.5); hold on;
grid on;
title('Gimbal vs. Time');
xlabel('Time [s]');
ylabel('Gimbal [rad]');

%% Integration
u_func = @(t) interp1(t_k(1:Nu), U_sol', t, "previous", "extrap")';
[t_cont, x_cont] = ode45(@(t, x) f(t, x, u_func(t), p_guess), t_k, x_0, tolerances);
x_cont = x_cont';

figure
plot(x_cont(1, :), x_cont(2, :));
grid on;
axis equal;
title('Integrated Intial Guess');
xlabel('X [m]');
ylabel('Y [m]');