import numpy as np
import cvxpy as cp
from cvxpygen import cpg
import time

A = np.random.rand(8, 173); # 8 x 173
y = np.random.rand(8, 1); # 8 x 1

x = cp.Variable((173, 1), name='x') # 173 x 1
A_param = cp.Parameter((8, 173), name = 'A')
y_param = cp.Parameter((8, 1), name = 'y') # 8 x 1
objective = cp.Minimize(cp.norm2(x))
constraints = [cp.norm2(A_param @ x - y_param) <= 0.1]

problem = cp.Problem(objective, constraints)

problem.param_dict["A"].value = A
problem.param_dict["y"].value = y

val = problem.solve(solver = "QOCO", verbose = True)


#cpg.generate_code(problem, code_dir=r'bigQOCO', solver = "QOCO")
# 
# 
#from bigQOCO.cpg_solver import cpg_solve
# 
#problem.register_solve('cpg', cpg_solve)
# 
# 
#t0 = time.time()
#val = problem.solve(method='cpg')
#t1 = time.time()
#print('\ncvxpy ecos_gen \nsolve time: %.3f ms with %.3f and %.5f ms solve' % (1000 * (t1 - t0), val, 1000 * problem.solution.attr["solve_time"]))
