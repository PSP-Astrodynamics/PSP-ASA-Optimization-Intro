function [X_sol, U_sol] = solving(w_tr_x, w_tr_u, w_vc, Delta_tol, x_tol, i_max, Delta, nx, nu, N, delta_t, x_guess, u_guess, glide_slope, T_max, gimbal_angle, m_dry, x_0, x_f, f, A, B, E, c, t_k, p_guess, tolerances)

x_ref = x_guess;
u_ref = u_guess;
p_ref = p_guess;

delta_x = 1000000000;

% PTR algorithm parameters
ptr_ops.iter_max = 25; % Maximum PTR iterations
ptr_ops.iter_min = 2; % Minimum PTR iterations
ptr_ops.Delta_min = 4e-5; % Convergence tolerance for sum total defects (difference between predicted and integrated final state) over each time interval
ptr_ops.w_vc = 1e4; % Virtual control penalty weight
ptr_ops.w_tr = 1e-1; % Trust region penalty weight
ptr_ops.delta_tol = 5e-3; % Convergence tolerance for change of new solution from reference solution
ptr_ops.alpha_x = 1; % State trust region multiplier
ptr_ops.alpha_u = 1e-6; % Control trust region multiplier

objective = @(x, u) sum(norms(u, 2, 1) * delta_t) / T_max;

[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_k, x_guess, u_guess, p_guess, tolerances);
total_defects(1) = sum(vecnorm(Delta, 2, 1));

i=1;

while (total_defects(i) > Delta_tol || delta_x > x_tol) && i <= i_max

    [X_sol, U_sol, V_sol] = solve_ptr_subproblem(objective, x_ref, u_ref, A_k, B_k, c_k, x_0, x_f, glide_slope, gimbal_angle, T_max, N, delta_t, ptr_ops);

    delta_x = sum(sum_square(X_sol - x_ref));

    % Linearize and discretize around new solution
    x_ref = X_sol;
    u_ref = U_sol;
    [A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_k, x_ref, u_ref, p_ref, tolerances);

    i = i+1;

    total_defects(i) = sum(vecnorm(Delta, 2, 1));

    fprintf("PTR Iter %.2g - Delta: %.3g, delta_x: %.3g\n", i, total_defects(i), delta_x);

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
    X_sol(1:(nx - 1), end) == x_f(1:(nx - 1)); % Terminal condition (final mass unconstrained)
    cvx_end
end