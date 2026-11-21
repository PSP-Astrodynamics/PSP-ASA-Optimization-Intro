function [X_sol, U_sol] = Solving_loop(w_tr_x, w_tr_u, w_vc, Delta_tol, x_tol, i_max, Delta, nx, nu, N, delta_t, x_guess, u_guess, glide_slope, T_max, gimbal_angle, m_dry, x_0, x_f, f, A, B, E, c, t_k, p_guess, tolerances)

x_ref = x_guess;
u_ref = u_guess;
p_ref = p_guess;

delta_x = 1000000000;

[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_k, x_guess, u_guess, p_guess, tolerances);
total_defects(1) = sum(vecnorm(Delta, 2, 1));

i=1;

while total_defects(i) > Delta_tol && delta_x > x_tol && i <= i_max

    cvx_begin quiet
        variable X_sol(nx, N)
        variable U_sol(nu, N - 1)
        variable V(nx, N-1)
        
        minimize( sum(norms(U_sol, 2, 1) * delta_t) + w_vc * sum(norms(V, 1, 1)) + w_tr_x * sum(sum_square(X_sol - x_guess)) + w_tr_u * sum(sum_square(U_sol - u_guess)))
        
        subject to
        % Dynamics constraints
        for k = 1 : (N - 1)
            X_sol(:, k + 1) == A_k(:, :, k) * X_sol(:, k) + B_k(:, :, k) * U_sol(:, k) + c_k(:, :, k) + V(:, k);
        end
        
        % Path constraints 
        norms(X_sol(1:2, :), 2, 1) <= X_sol(2, :)/cos(glide_slope); % Max distance constraint
        norms(U_sol, 2, 1) <= T_max; % Max control magnitude constraint
        norms(U_sol, 2, 1) <= U_sol(1, :)*cos(gimbal_angle) % Max gimbal
        X_sol(7, :) >= m_dry %
        
        % State Boundary Conditions
        X_sol(:, 1) == x_0; % Initial condition
        X_sol(:, end) == x_f; % Terminal condition
    cvx_end

    delta_x = sum(sum_square(X_sol - x_ref));

    % Linearize and discretize around new solution
    x_ref = X_sol;
    u_ref = U_sol;
    [A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, E, c, N, t_k, x_ref, u_ref, p_ref, tolerances);

    i = i+1;

    total_defects(i) = sum(vecnorm(Delta, 2, 1));

    fprintf("PTR Iter %.2g - Delta: %.3g, delta_x: %.3g\n", i, total_defects(i), delta_x);

end