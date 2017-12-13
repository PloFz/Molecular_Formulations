import general_functions as gf, numpy as np, bempp.api, inspect, time, PBL, os, sys

def solvation_energy(mol_name, mesh_density, ep_in=4., ep_ex=80., kappa=0.125, solver_tol=1e-3,
                     formulation='stern_d', stern_radius=1.4, info=False, phi_info=False):
    '''
    mol_name  : molecule name to call .pqr & .msh files from its directories
    mesh_density : mesh density from msms cration
    ep_in : interior electrostatic permittivity
    ep_ex : exterior electrostatic permittivity
    kappa : kappa constant for Poisson-Boltzmann field
    '''
    total_time = time.time()

    mesh_directory = 'Molecule/' + mol_name + '/mesh'
    mesh_file_in = '{}/{}_d{:04.1f}.msh'.format(mesh_directory, mol_name, mesh_density)
    mesh_file_ex = '{}/{}_d{:04.1f}_strn-pr{:04.1f}.msh'.format(mesh_directory, mol_name, mesh_density, stern_radius)

    grid_in = bempp.api.import_grid(mesh_file_in)
    grid_ex = bempp.api.import_grid(mesh_file_ex)

    # Define spaces for both boundary phi & dphi
    dirichl_space_in = bempp.api.function_space(grid_in, "DP", 0)
    neumann_space_in = bempp.api.function_space(grid_in, "DP", 0)
    dirichl_space_ex = bempp.api.function_space(grid_ex, "DP", 0)
    neumann_space_ex = bempp.api.function_space(grid_ex, "DP", 0)

    q, x_q = gf.read_pqr(mol_name)

    print "\nNew Calculation"
    print "Molecule name : " + mol_name
    print "Formulation   : " + formulation
    #print "Mesh file: {}".format(mesh_file)
    print "Mesh density : {:5.2f}".format(mesh_density)
    print "Stern Radius : {:5.2f}".format(stern_radius)
    print "Number of elements in  : {0}".format(dirichl_space_in.global_dof_count)
    print "Number of elements out : {0}".format(dirichl_space_ex.global_dof_count)

    matrix_time = time.time()
    if formulation == 'stern_d':
        A, rhs = stern_formulation(dirichl_space_in, neumann_space_in, 
                                   dirichl_space_ex, neumann_space_ex, 
                                   ep_in, ep_ex, q, x_q, kappa)
    elif formulation == 'asc':
        A, rhs = stern_asc(neumann_space_in, dirichl_space_ex, neumann_space_ex, 
                           ep_in, ep_ex, q, x_q, kappa)
    matrix_time = time.time() - matrix_time

    # Preconditioning
    #A_prec = gf.inverse_block_diagonals(A)

    # Solver GMRES
    print "Solving system...\n"
    solver_time = time.time()
    from scipy.sparse.linalg import gmres
    global array_it, array_frame, it_count
    array_it, array_frame, it_count = np.array([]), np.array([]), 0
    x, _ = gmres(A, rhs, M=A_prec, callback=iteration_counter, tol=solver_tol, maxiter=1000, restart = 1000)
    solver_time = time.time() - solver_time
    print "\nSolved in {} iterations".format(it_count)
    print "Residual tolerance : {0:1.2e}".format(solver_tol)

    if phi_info:
        phi_file = open('{}_{}_{:5.2f}'.format(mol_name, formulation, mesh_density),'w')
        if formulation=='asc':
            np.savetxt(phi_file, x)
        else:
            np.savetxt(phi_file, x[dirichl_space_in.global_dof_count:])

    from bempp.api.operators import potential
    if formulation == 'stern_d':
        # The following grid function stores the computed boundary data of the total field.
        p1, p2 = np.split(x[:(dirichl_space_in.global_dof_count + neumann_space_in.global_dof_count)], 2) 
        total_dirichl_in = bempp.api.GridFunction(dirichl_space_in, coefficients=p1)
        total_neumann_in = bempp.api.GridFunction(neumann_space_in, coefficients=p2)

        # Calculate potentials in the coordinates of the atoms
        slp_ev = potential.laplace.single_layer(neumann_space_in, x_q.transpose())
        dlp_ev = potential.laplace.double_layer(dirichl_space_in, x_q.transpose())

        # Evaluate potential at charges position & total dissolution energy
        phi_q = slp_ev*total_neumann_in - dlp_ev*total_dirichl_in

    elif formulation == 'asc':
        sigma_in = np.array((ep_in/ep_ex - 1.)*x[:neumann_space_in.global_dof_count])
        total_sigma_in = bempp.api.GridFunction(neumann_space_in, coefficients=sigma_in)
	slp_ev = potential.laplace.single_layer(neumann_space_in, x_q.transpose())
        phi_q =  slp_ev*total_sigma_in

    total_energy = 2*np.pi*332.064*np.sum(q*phi_q).real

    total_time = time.time() - total_time

    print "\nAssamble time : {:5.2f}".format(matrix_time)
    print "Solver time   : {:5.2f}".format(solver_time)
    print "Total time    : {:5.2f}".format(total_time)
    print "\nTotal dissolution Energy: {:8.3f}".format(total_energy)

    if info:
        info_dict = {}
        info_dict['mol_name'] = mol_name
        info_dict['mesh_density'] = mesh_density
        info_dict['n_of_elements'] = dirichl_space_in.global_dof_count
        info_dict['n_of_elements_ex'] = dirichl_space_ex.global_dof_count
        info_dict['formulation'] = formulation
        info_dict['energy'] = total_energy
        info_dict['assemble_time'] = matrix_time
        info_dict['solver_time'] = solver_time
        info_dict['total_time'] = total_time
        info_dict['iterations'] = it_count
        info_dict['stern_radius'] = stern_radius

        gf.save_log(mol_name, info_dict)

    return total_energy

def stern_formulation(dirichl_space_in, neumann_space_in, dirichl_space_ex, neumann_space_ex, 
                      ep_in, ep_ex, q, x_q, kappa):

    # Functions to proyect the carges potential to the boundary with constants
    def green_func(x, n, domain_index, result):
        result[:] = np.sum(q/np.linalg.norm( x - x_q, axis=1 ))/(4.*np.pi*ep_in)

    print "\nProjecting charges over surface..."
    charged_grid_fun  = bempp.api.GridFunction(dirichl_space_in, fun=green_func)

    rhs = np.concatenate([charged_grid_fun.coefficients, 
                          np.zeros(neumann_space_in.global_dof_count),
                          np.zeros(dirichl_space_ex.global_dof_count), 
                          np.zeros(neumann_space_ex.global_dof_count)])

    print "Defining operators..."
    # OPERATOR FOR INTERNAL SURFACE
    from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz
    idn_in  = sparse.identity(dirichl_space_in, dirichl_space_in, dirichl_space_in)

    # Ec 1
    slp_in = laplace.single_layer(neumann_space_in, dirichl_space_in, dirichl_space_in)
    dlp_in = laplace.double_layer(dirichl_space_in, dirichl_space_in, dirichl_space_in)
    # Ec 2
    # adj_1T1 = laplace.single_layer(neumann_space_in, dirichl_space_in, dirichl_space_in)
    # adj_1T1 = laplace.double_layer(dirichl_space_in, dirichl_space_in, dirichl_space_in)

    slp_2T1 = laplace.single_layer(neumann_space_ex, dirichl_space_in, dirichl_space_in)
    dlp_2T1 = laplace.double_layer(dirichl_space_ex, dirichl_space_in, dirichl_space_in)

    # OPERATOR FOR EXTERNAL SURFACE
    idn_ex = sparse.identity(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex)
    # Internal Boudary
    slp_1T2 = laplace.single_layer(neumann_space_in, dirichl_space_ex, dirichl_space_ex)
    # dlp_1T2 = laplace.double_layer(dirichl_space_in, dirichl_space_ex, dirichl_space_ex)

    slp_2T2 = laplace.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex)
    dlp_2T2 = laplace.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex)
    # External Boundary
    slp_ex = modified_helmholtz.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)
    dlp_ex = modified_helmholtz.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)

    ep = (ep_in/ep_ex)

    print "Creating operators..."
    # Matrix Assemble
    blocked = bempp.api.BlockedOperator(4, 4)
    blocked[0, 0] = .5*idn_in + dlp_in
    blocked[0, 1] = -slp_in
    # blocked[0, 2] = 0
    # blocked[0, 3] = 0

    # Original formulation
    blocked[1, 0] = .5*idn_in - dlp_in   # dlp_in = dlp_1T1
    blocked[1, 1] =  ep*slp_in           # slp_in = slp_1T1
    blocked[1, 2] =  dlp_2T1
    blocked[1, 3] = -slp_2T1

    # blocked[2, 0] = -dlp_1T2    ## eliminar**
    blocked[2, 1] =  ep*slp_1T2
    blocked[2, 2] = .5*idn_ex + dlp_2T2
    blocked[2, 3] = -slp_2T2

    # blocked[3, 0] = 0
    # blocked[3, 1] = 0
    blocked[3, 2] = .5*idn_ex - dlp_ex
    blocked[3, 3] = slp_ex
    A = blocked.strong_form()

    return A, rhs


def stern_asc(sigma_space_in, dirichl_space_ex, neumann_space_ex, 
              ep_in, ep_ex, q, x_q, kappa):

    # Functions to proyect the carges potential to the boundary with constants
    def d_green_func(x, n, domain_index, result):
        const = -1./(4.*np.pi*ep_in)
        result[:] = const*np.sum(q*np.dot( x - x_q, n )/(np.linalg.norm( x - x_q, axis=1 )**3))

    print "\nProjecting charges over surface..."
    charged_grid_fun  = bempp.api.GridFunction(sigma_space_in, fun=d_green_func)

    rhs = np.concatenate([charged_grid_fun.coefficients,
                          np.zeros(dirichl_space_ex.global_dof_count),
                          np.zeros(neumann_space_ex.global_dof_count)])

    print "Defining operators..."
    # OPERATORS FOR INTERNAL SURFACE
    from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz
    idn_in  = sparse.identity(sigma_space_in, sigma_space_in, sigma_space_in)

    adj_1T1 = laplace.adjoint_double_layer(sigma_space_in, sigma_space_in, sigma_space_in)
    hyp_2T1 = laplace.hypersingular(dirichl_space_ex, sigma_space_in, sigma_space_in)
    adj_2T1 = laplace.adjoint_double_layer(neumann_space_ex, sigma_space_in, sigma_space_in)

    # OPERATORS FOR EXTERNAL SURFACE
    idn_ex = sparse.identity(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex)

    slp_1T2 = laplace.single_layer(sigma_space_in, dirichl_space_ex, dirichl_space_ex)
    dlp_2T2 = laplace.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex)
    slp_2T2 = laplace.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex)

    dlp_ex = modified_helmholtz.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)
    slp_ex = modified_helmholtz.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)

    ep = (ep_in/ep_ex)

    print "Creating operators..."
    # Matrix Assemble
    blocked = bempp.api.BlockedOperator(3, 3)
    blocked[0, 0] = idn_in + (ep - 1.)*adj_1T1
    blocked[0, 1] = -hyp_2T1
    blocked[0, 2] = -adj_2T1

    blocked[1, 0] = ep*slp_1T2
    blocked[1, 1] = .5*idn_ex + dlp_2T2
    blocked[1, 2] = -slp_2T2

    #blocked[2, 0] = 0
    blocked[2, 1] = .5*idn_ex - dlp_ex
    blocked[2, 2] = slp_ex

    A = blocked.strong_form()

    return A, rhs

array_it, array_frame, it_count = np.array([]), np.array([]), 0
def iteration_counter(x):
    global array_it, array_frame, it_count
    it_count += 1
    frame, array_it = inspect.currentframe().f_back, np.append(array_it, it_count)
    array_frame = np.append(array_frame, frame.f_locals["resid"])
    sys.stdout.flush()
    sys.stdout.write("Iter: {0} - Error: {1:.2E}    \r".format(it_count, frame.f_locals["resid"]))
