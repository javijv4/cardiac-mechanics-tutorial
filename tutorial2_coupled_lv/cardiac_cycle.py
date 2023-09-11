#!/usr/bin/env python3
"""
Cardiac Mechanics Tutorial w/ AMBIT
Cardiac Biomechanics Lab @ UofM
Prof. David Nordsletten

@author: Javiera Jilberto Vallejos (jilberto@umich.edu)
"""
import os
import ambit_fe
import numpy as np


def main():
    data_path = 'data'
    out_path = 'out_py'

    if not os.path.exists(out_path): os.mkdir(out_path)     # Make sure the output path exists

    IO_PARAMS            = {'problem_type'          : 'solid_flow0d',            # Defines the type of problem to be solved
                            # Inputs
                            'mesh_domain'           : data_path + '/mesh.xdmf',   # Path to mesh
                            'mesh_boundary'         : data_path + '/mt.xdmf',     # Path to boundary file
                            'order_fib_input'       : 1,                         # Define the interpolation order of the fibers
                            'fiber_data'            : [data_path + '/fiber.txt', data_path + '/sheet.txt', data_path + '/sheetnormal.txt'],  # Path to fibers (f,s,n)
                            # Output options
                            'output_path'           : out_path,           # Path where the results will be saved
                            'write_results_every'   : 1,                 # Saving steps every n time steps
                            'results_to_write'      : ['displacement', 'pressure'],  # Which results to write
                            'simname'               : 'test'}             # Name of the simulation, all the results will have this name

    SOLVER_PARAMS        = {'solve_type'            : 'direct',          # Linear Algebra solver (direct or iterative)
                            'tol_res'               : [1.0e-8,1.0e-8,1.0e-6],   # Residual tolerance for [displacements, pressure, 0D problem]
                            'tol_inc'               : [1.0e-8,1.0e-8,1.0e-6],   # Increment tolerance for [displacements, pressure, 0D problem]
                            }

    TIME_PARAMS_SOLID    = {'maxtime'               : 10.0,      # Final time. In this case is 10 sec.
                            'numstep'               : 5000,      # Solving the 10 sec using 5000 timesteps, i.e., dt=0.002
                            'numstep_stop'          : 5000,      # If you want to stop the simulation before set this to whatever timestep you want to stop
                            'timint'                : 'ost',     # Time integration scheme. This is the trapezoidal rule
                            'theta_ost'             : 1.0}       # Trapezoidal rule parameter. 1.0 means backward euler.

    TIME_PARAMS_FLOW0D   = {'timint'                : 'ost',      # Use trapezoidal rule
                            'theta_ost'             : 1.0,        # Trapezoidal rule parameter. 1.0 means backward euler.
                            'initial_conditions'    : init()}     # Initial conditions of the 0D problem

    MODEL_PARAMS_FLOW0D  = {'modeltype'             : 'syspul',
                            'parameters'            : param(),
                            'chamber_models'        : {'lv' : {'type' : '3D_solid'},
                                                       'rv' : {'type' : '0D_elast', 'activation_curve' : 1},
                                                       'la' : {'type' : '0D_elast', 'activation_curve' : 2},
                                                       'ra' : {'type' : '0D_elast', 'activation_curve' : 2}}}

    FEM_PARAMS           = {'order_disp'            : 1,
                            'order_pres'            : 1,
                            'quad_degree'           : 4,       # Quadrature degree
                            'incompressible_2field' : True}    # If using or not the pressure field

    COUPLING_PARAMS      = {'surface_ids'           : [[3]],
                            'surface_p_ids'         : [[3]],
                            'coupling_quantity'     : ['flux'],             # The change in volume (flux) will be the coupling quantity
                            'coupling_type'         : 'monolithic_direct'}  # How is the coupled system going to be solved.

    MATERIALS            = {'MAT1' : {'holzapfelogden_dev' : {'a_0' : 0.4, 'b_0' : 3.2, 'a_f' : 1.0, 'b_f' : 5, 'a_s' : 0., 'b_s' : 0.1, 'a_fs' : 0., 'b_fs' : 0.1, 'fiber_comp' : False},
                                      'active_fiber'      : {'sigma0' : 60.0, 'alpha_max' : 15.0, 'alpha_min' : -20.0, 'activation_curve' : 1, 'frankstarling' : False, 'amp_min' : 1., 'amp_max' : 1.7, 'lam_threslo' : 1.01, 'lam_maxlo' : 1.15, 'lam_threshi' : 999., 'lam_maxhi' : 9999.},
                                      'inertia'           : {'rho0' : 1.0e-6},
                                      'visco_green'       : {'eta' : 0.1}}}



    # Read data for activation and define a function to evaluate it.
    activation = np.loadtxt(data_path + '/activation.txt')

    from scipy.interpolate import interp1d
    func_act = interp1d(activation[:,0], activation[:,1])

    # Define atrial and ventricle activation curves
    class time_curves():
        def tc1(self, t): # ventricle activation
            K = 5.
            t_contr, t_relax = 0.139, 0.535
            alpha_max = MATERIALS['MAT1']['active_fiber']['alpha_max']
            alpha_min = MATERIALS['MAT1']['active_fiber']['alpha_min']
            c1 = t_contr + alpha_max/(K*(alpha_max-alpha_min))
            c2 = t_relax - alpha_max/(K*(alpha_max-alpha_min))

            # Diss Hirschvogel eq. 2.101
            return (K*(t-c1)+1.)*((K*(t-c1)+1.)>0.) - K*(t-c1)*((K*(t-c1))>0.) - K*(t-c2)*((K*(t-c2))>0.) + (K*(t-c2)-1.)*((K*(t-c2)-1.)>0.)


        def tc2(self, t): # atrial activation

            act_dur = 2.*param()['t_ed']
            t0 = 0.

            if t >= t0 and t <= t0 + act_dur:
                return 0.5*(1.-np.cos(2.*np.pi*(t-t0)/act_dur))
            else:
                return 0.0


    BC_DICT              = { 'dirichlet' : [{'id' : [2], 'dir' : 'all', 'val' : 0.}] }

    # problem setup
    problem = ambit_fe.ambit_main.Ambit(IO_PARAMS, [TIME_PARAMS_SOLID, TIME_PARAMS_FLOW0D], SOLVER_PARAMS, FEM_PARAMS, [MATERIALS, MODEL_PARAMS_FLOW0D], BC_DICT, time_curves=time_curves(), coupling_params=COUPLING_PARAMS)

    # problem solve
    problem.solve_problem()

    return




# syspulcap circulation model initial condition and parameter dicts...

def init():
    return {'p_at_l_0' : 1.213325556608718,
            'p_ar_sys_0' : 9.15625635011591,
            'p_ven_sys_0' : 1.8954654756677949,
            'p_at_r_0' : 0.5026616990281122,
            'p_v_r_0' : 0.4470328846734767,
            'p_ar_pul_0' : 1.735368268177882,
            'p_ven_pul_0' : 1.5577155509540659,
            'q_vin_l_0' : 0.,
            'q_vout_l_0' : 0.,
            'q_ar_sys_0' : 0.,
            'q_ar_pul_0' : 0.,
            'q_ven_sys_0' : 0.,
            'q_ven_pul_0' : 0.,
            'q_vin_r_0' : 0.,
            'q_vout_r_0' : 0.,
            'p_v_l_0' : 0.,
            }



def param():
    return {'R_ar_sys' : 0.00012173913043478261,
            'C_ar_sys' : 12321.42857142857,
            'L_ar_sys' : 6.67e-07,
            'Z_ar_sys' : 6.086956521739131e-06,
            'R_ar_pul' : 7.699999999999999e-06,
            'C_ar_pul' : 19714.285714285714,
            'L_ar_pul' : 1e-08,
            'Z_ar_pul' : 0.,
            'R_ven_sys' : 1.8260869565217376e-05,
            'C_ven_sys' : 369642.8571428571,
            'L_ven_sys' : 1e-08,
            'R_ven_pul' : 6.299999999999999e-06,
            'C_ven_pul' : 49285.71428571428,
            'L_ven_pul' : 1e-08,
            # atrial elastances
            'E_at_max_l' : 39e-06,
            'E_at_min_l' : 1.8e-05,
            'E_at_max_r' : 12e-06,
            'E_at_min_r' : 8e-06,
            # ventricular elastances
            'E_v_max_l' : 20e-05,
            'E_v_min_l' : 1e-05,
            'E_v_max_r' : 5e-05,
            'E_v_min_r' : 4e-06,
            # valve resistances
            'R_vin_l_min' : 1.0e-6,
            'R_vin_l_max' : 1.0e1,
            'R_vout_l_min' : 1.0e-6,
            'R_vout_l_max' : 1.0e1,
            'R_vin_r_min' : 1.0e-6,
            'R_vin_r_max' : 1.0e1,
            'R_vout_r_min' : 1.0e-6,
            'R_vout_r_max' : 1.0e1,
            # timings
            't_ed' : 0.139,
            't_es' : 0.535,
            'T_cycl' : 1.0,
            }




if __name__ == "__main__":

    main()
