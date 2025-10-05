import numpy as np
import cvxpy as cp
from cvxpygen import cpg
import time
import pickle

save = False;
load = True;

if not load:
    # Input from matlab
    x_ref = np.array(x_ref)
    u_ref = np.array(u_ref)

    x_0 = np.array(x_0).reshape(-1, 1)
    x_f = np.array(x_f).reshape(-1, 1)

    # Problem parameters
    # T_min, T_max, tau_max, alpha, glideslope_max_angle, gimbal_max_angle, pitch_max, angvel_max, time_min_max_thrust, max_gimbal_rate
    T_min = np.array(params[0])
    T_max = np.array(params[1])
    glideslope_max_angle = np.array(params[2])
    gimbal_max_angle = np.array(params[3])
    S_x = np.reshape(np.array(params[0 : 6]), [6, 1]) 
    S_u = np.reshape(np.array(params[6 : 9]), [3, 1])
    c_x = np.reshape(np.array(params[9 : 15]), [6, 1])
    c_u = np.reshape(np.array(params[15 : 18]), [3, 1])
    N = int(N);
    delta_t = np.array(delta_t)

    # PTR parameters
    w_vc = np.array(w_vc)
    w_tr = np.array(w_tr)

    ### Discretization
    A_k = np.array(A_k)
    B_k_minus = np.array(B_k_minus)
    B_k_plus = np.array(B_k_plus)
    c_k = np.array(c_k)
else:
    loaded_data = np.load(r'Planar Rocket Landing\CVXPyGEN Example\Deterministic_3DoF_with_mass_constant_input_FOH.npz', allow_pickle=True)
    x_ref = loaded_data["x_ref"]
    u_ref = loaded_data["u_ref"]

    x_0 = loaded_data["x_0"]
    x_f = loaded_data["x_f"]

    # Problem parameters
    # T_min, T_max, glideslope_max_angle, gimbal_max_angle
    T_min = loaded_data["T_min"]
    T_max = loaded_data["T_max"]
    glideslope_max_angle = loaded_data["glideslope_max_angle"]
    gimbal_max_angle = loaded_data["gimbal_max_angle"]
    N = loaded_data["N"]
    delta_t = loaded_data["delta_t"]
    S_x = loaded_data["S_x"]
    S_u = loaded_data["S_u"]
    c_x = loaded_data["c_x"]
    c_u = loaded_data["c_u"]

    # PTR parameters
    w_vc = loaded_data["w_vc"]
    w_tr = loaded_data["w_tr"]

    # Discretization
    A_k = loaded_data["A_k"]
    B_k_minus = loaded_data["B_k_minus"]
    B_k_plus = loaded_data["B_k_plus"]
    c_k = loaded_data["c_k"]

    # Saved Problem
    dbfile = open(r'Planar Rocket Landing\CVXPyGEN Example\Deterministic_3DoF_with_mass_constant_prob_FOH', 'rb')    
    problem = pickle.load(dbfile)

    dbfile.close()

t_k = np.linspace(0, N * delta_t, N + 1)

nx = 6
nu = 3

def virtual_control_cost_func(w_vc, V, v_0, v_N):
    J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.norm(v_0, 1) + cp.norm(v_N, 1))
    return J_vc

def trust_region_cost_func(w_tr, eta):
    J_tr = w_tr @ eta.T
    return J_tr

def trust_region_cost_func_explicit(w_tr, x, x_ref, u, u_ref):
    J_tr = (cp.sum(cp.square(x - x_ref), 0) + cp.sum(cp.square(u - u_ref), 0)) @ w_tr.T
    return J_tr

def virtual_control_cost_func_np(w_vc, V, v_0, v_N):
    J_vc = w_vc * (np.sum(np.linalg.norm(V, 1, 0)) + np.linalg.norm(v_0, 1) + np.linalg.norm(v_N, 1))
    return J_vc

def scale_x(x, S_x, c_x):
    x_scaled = (x - c_x) / S_x
    return x_scaled 


## Create Convex Problem
# Define variables
if True:#str(type(problem)) == "<class 'array.array'>" or str(type(problem)) == "<class 'numpy.ndarray'>":
    #print("Creating problem")
    X = cp.Variable((nx, N + 1), name='X')
    U = cp.Variable((nu, N + 1), name='U')
    #eta = cp.Variable((1, N + 1), name = 'eta')
    V = cp.Variable((nx, N + 0), name = 'V')
    v_0 = cp.Variable((nx, 1), name = 'v_0')
    v_N = cp.Variable((nx, 1), name = 'v_N')

    # Define parameters
    ck_param = cp.Parameter((nx, N + 0), name = 'ck')
    x_0_param = cp.Parameter((nx, 1), name = 'x_0')
    x_f_param = cp.Parameter((nx, 1), name = 'x_f')

    x_ref_param = cp.Parameter((nx, N + 1), name = 'x_ref')
    u_ref_param = cp.Parameter((nu, N + 1), name = 'u_ref')

    # Create unscaled variables
    X_unscaled = cp.multiply(S_x, X) + c_x;
    U_unscaled = cp.multiply(S_u, U) + c_u;

    # Define objective
    objective = cp.Minimize(cp.sum(U_unscaled[2, 0:-1] + U_unscaled[2, 1:]) / 2 * delta_t) # Min total thrust magnitude (min fuel)
    virtual_control_cost = cp.Minimize(virtual_control_cost_func(w_vc, V, v_0, v_N))
    #trust_region_cost = cp.Minimize(trust_region_cost_func(w_tr, eta)) # 
    trust_region_cost = cp.Minimize(trust_region_cost_func_explicit(w_tr, X, x_ref_param, U, u_ref_param))
    augmented_objective = objective + virtual_control_cost + trust_region_cost

    # Define constraints
    Ak_params = []
    Bk_minus_params = []
    Bk_plus_params = []
    dynamics_constraints = []
    for k in range(N):
        Ak_params.append(cp.Parameter((nx, nx), name = "Ak_" + str(k)))
        Bk_plus_params.append(cp.Parameter((nx, nu), name = "Bk_plus_" + str(k)))
        Bk_minus_params.append(cp.Parameter((nx, nu), name = "Bk_minus_" + str(k)))

        dynamics_constraints.append(X[:, k + 1] == scale_x(Ak_params[k] @ X_unscaled[:, k] + Bk_minus_params[k] @ U_unscaled[:, k] + Bk_plus_params[k] @ U_unscaled[:, k + 1] + ck_param[:, k] + V[:, k], S_x, c_x))

    # Min thrust, max thrust, glideslope, max gimbal, max roll torque, max angular velocity, thrust rate
    constraints = dynamics_constraints + [ # Dynamics constraint
                   T_min <= U_unscaled[2, :], # Thrust min constraint
                   U_unscaled[2, :] <= T_max, # Thrust max constraint
                   cp.norm2(X_unscaled[0:2, :], 0) <= X_unscaled[1, :] / np.cos(glideslope_max_angle), # Glideslope constraint
                   U_unscaled[2, :] <= U_unscaled[0, :] / np.cos(gimbal_max_angle),
                   cp.norm2(U_unscaled[0:2, :], 0) <= U_unscaled[2, :], # Lcvx constraint
                   X_unscaled[:, 0] + v_0 == x_0_param.flatten(), # Initial condition constraint
                   X_unscaled[:, -1] + v_N == x_f_param.flatten()] # Terminal condition constraint
                   #cp.sum(cp.sum_squares(X - x_ref_param)) + cp.sum(cp.sum_squares(U - u_ref_param)) <= eta]

    # Define problem
    problem = cp.Problem(augmented_objective, constraints)


# Save input data so it can be debugged in python
if save:
    np.savez(r'Planar Rocket Landing\CVXPyGEN Example\Deterministic_3DoF_with_mass_constant_input_FOH.npz', A_k = A_k, B_k_minus = B_k_minus, B_k_plus = B_k_plus, c_k = c_k, w_vc = w_vc, w_tr = w_tr, T_min = T_min, T_max = T_max, glideslope_max_angle = glideslope_max_angle, gimbal_max_angle = gimbal_max_angle, delta_t = delta_t, N = N, x_0 = x_0, x_f = x_f, x_ref = x_ref, u_ref = u_ref, S_x = S_x, S_u = S_u, c_x = c_x, c_u = c_u)
    
    dbfile = open(r'Planar Rocket Landing\CVXPyGEN Example\Deterministic_3DoF_with_mass_constant_prob_FOH', 'ab')
    
    # source, destination
    pickle.dump(problem, dbfile)                    
    dbfile.close()

# Set parameters
for k in range(N):
    problem.param_dict["Ak_" + str(k)].value = A_k[:, :, k]
    problem.param_dict["Bk_minus_" + str(k)].value = B_k_minus[:, :, k]
    problem.param_dict["Bk_plus_" + str(k)].value = B_k_plus[:, :, k]

problem.param_dict["ck"].value = c_k[:, 0, 0:N]

problem.param_dict["x_0"].value = x_0
problem.param_dict["x_f"].value = x_f

problem.param_dict["x_ref"].value = x_ref
problem.param_dict["u_ref"].value = u_ref

# Solve
t0 = time.time()
val = problem.solve(solver = "CLARABEL", verbose = True)#solver = "ECOS")
t1 = time.time()
#print('\nCVXPY Clarabel\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))
val = problem.solve(solver = "QOCO", verbose = True)#solver = "ECOS")

# cpg.generate_code(problem, code_dir=r'Planar Rocket Landing\CVXPyGEN Example\Deterministic_3DoF_with_mass_constant_FOH_QOCO', solver = "QOCO")
# 
# 
# from Deterministic_3DoF_with_mass_constant_FOH_QOCO.cpg_solver import cpg_solve
# 
# problem.register_solve('cpg', cpg_solve)
# 
# 
# t0 = time.time()
# val = problem.solve(method='cpg')
# t1 = time.time()
# print('\ncvxpy ecos_gen \nsolve time: %.3f ms with %.3f and %.5f ms solve' % (1000 * (t1 - t0), val, 1000 * problem.solution.attr["solve_time"]))

# Extract Solution
X_sol = problem.var_dict['X'].value
U_sol = problem.var_dict['U'].value

eta = np.sum(np.square(X_sol - x_ref), 0) + np.sum(np.square(U_sol - u_ref), 0) #problem.var_dict['eta'].value
V = problem.var_dict['V'].value
v_0 = problem.var_dict['v_0'].value
v_N = problem.var_dict['v_N'].value

solve_status = problem.status