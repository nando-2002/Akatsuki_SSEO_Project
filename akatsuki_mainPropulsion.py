import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numba import njit, objmode, prange
import matplotlib.pyplot as plt
from CEA_Wrap import Fuel, Oxidizer, RocketProblem # pyright: ignore[reportMissingImports]

# Akatsuki Main Propulsion System 
# 500 N Thrust 
# 0.8 MPa
# N2H4 and NTO


# Propellants
hydrazine = Fuel("N2H4(L)")
nto = Oxidizer("N2O4(L)")

problem_1 = RocketProblem(pressure = 6.9, materials = [hydrazine, nto], o_f = 0.8, sup=50)
problem_1.set_pressure_units("bar")
results = problem_1.run()

of_range = np.linspace(0.1, 2, 100)
area_range = np.linspace(1, 100, 100) 
results_arr = np.zeros((len(of_range), len(area_range)))

@njit(parallel=True)
def compute_results_objmode(of_range, area_range, results_arr):
    for i in prange(len(of_range)):
        for j in prange(len(area_range)):
            # Initialize a variable to hold the result
            isp_val = 0.0 
            
            # Entering the Python bridge
            # 'isp_val' must match the variable name inside the block
            with objmode(isp_val='float64'):
                # Regular Python code goes here
                # Numba cannot see this, it just waits for it to finish
                prob = RocketProblem(
                    pressure=6.9, 
                    materials=[hydrazine, nto], 
                    o_f=of_range[i], 
                    sup=area_range[j]
                )
                res = prob.run()
                isp_val = res.isp
            
            # Back in Numba mode
            results_arr[i, j] = isp_val

compute_results_objmode(of_range, area_range, results_arr)

""" for i in range(len(of_range)): 
    for j in range(len(area_range)): 
        curr_problem = RocketProblem(
            pressure = 6.9, 
            materials = [hydrazine, nto], 
            o_f = of_range[i], 
            sup=area_range[j])
        exit_isp = curr_problem.run()

        results_arr[i, j] = exit_isp.isp
 """
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(area_range, of_range)
ax.plot_surface(X, Y, results_arr, cmap='viridis')

ax.set_xlabel('Supersonic Area Ratio Range')
ax.set_ylabel('O/F Range')
ax.set_zlabel('ISP')
plt.savefig("hello2")
""" 
print("Akatsuki Main Thruster")
exit_pressure    = results.p
chamber_pressure = results.c_p
exit_cp          = results.cp
chamber_cp       = results.c_cp
exit_isp         = results.isp
throat_isp       = results.t_isp
print("Exit Pressure (bar):", exit_pressure," | Chamber Pressure: ", chamber_pressure)
print("Exit Cp (kJ/(kg*K):", exit_cp," | Chamber Cp: ", chamber_cp)
print("Exit Isp (s):", exit_isp," | Throat Isp: ", throat_isp)
"""