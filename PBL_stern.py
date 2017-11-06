import numpy as np
import bempp.api
import time
import PBL
import os

from scipy.sparse.linalg import gmres

def PBL_stern(mol_name, mesh_density, ep_in=4., ep_ex=80., kappa=0.125, 
                formulation='stern_d', probe_radius=1.0, info=False):
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
    mesh_file_ex = '{}/{}_d{:04.1f}_strn-pr{:04.1f}.msh'.format(mesh_directory, mol_name, mesh_density)

    grid_in = bempp.api.import_grid(mesh_file_in)
    grid_ex = bempp.api.import_grid(mesh_file_ex)

    # Define spaces for both boundary phi & dphi
    dirichl_space_in = bempp.api.function_space(grid_in, "DP", 0)
    neumann_space_in = bempp.api.function_space(grid_in, "DP", 0)
    dirichl_space_ex = bempp.api.function_space(grid_ex, "DP", 0)
    neumann_space_ex = bempp.api.function_space(grid_ex, "DP", 0)

    # number of elemnts in both grids
    n_elmts_in = dirichl_space_in.global_dof_count
    n_elmts_ex = dirichl_space_ex.global_dof_count

    q, x_q = PBL.read_pqr(mol_name)

    def charges_fun(x, n, domain_index, result):
        result[:] = np.sum(q/np.linalg.norm( x - x_q, axis=1 ))/(4.*np.pi*ep_in)

    print "Proyecting charges over surface..."
    charged_grid_fun  = bempp.api.GridFunction(dirichl_space_in, fun=charges_fun)

    rhs = np.concatenate([charged_grid_fun.coefficients,
			              np.zeros(n_elmts_in),
			              np.zeros(n_elmts_ex),
			              np.zeros(n_elmts_ex)])

    print "Assembling Matrix"
    matrix_time = time.time()

    # Define Operators
    from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

    # OPERATOR FOR INTERNAL SURFACE
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
    dlp_1T2 = laplace.double_layer(dirichl_space_in, dirichl_space_ex, dirichl_space_ex)

    slp_2T2 = laplace.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex)
    dlp_2T2 = laplace.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex)
    # External Boundary
    slp_ex = modified_helmholtz.single_layer(neumann_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)
    dlp_ex = modified_helmholtz.double_layer(dirichl_space_ex, dirichl_space_ex, dirichl_space_ex, kappa)

    ep = (ep_ex/ep_in)
    # Matrix Assembly
    blocked = bempp.api.BlockedOperator(4, 4)
    blocked[0, 0] = 0.5*ident_in + dlp_in
    blocked[0, 1] = - ep*slp_in
    #blocked[0, 2] = 0
    #blocked[0, 3] = 0

    #Original formulation
    blocked[1, 0] = 0.5*ident_in - dlp_in   # dlp_in = dlp_1T1
    blocked[1, 1] =  slp_in                 # slp_in = slp_1T1
    blocked[1, 2] =  dlp_2T1
    blocked[1, 3] = -slp_2T1

    blocked[2, 0] = -dlp_1T2
    blocked[2, 1] =  slp_1T2
    blocked[2, 2] = 0.5*ident_ex + dlp_2T2
    blocked[2, 3] = -slp_2T2

    #blocked[3, 0] = 0
    #blocked[3, 1] = 0
    blocked[3, 2] = 0.5*ident_ex - dlp_ex
    blocked[3, 3] = slp_ex
    A = blocked.strong_form()

    matrix_time = time.time() - matrix_time
    print "Assamble time: {:5.2f}".format(matrix_time)

    # Solver GMRES
    solver_time = time.time()
    global array_it, array_frame, it_count
    array_it, array_frame, it_count = np.array([]), np.array([]), 0
    x, info = gmres(A, rhs, callback=PBL.iteration_counter, tol=1e-3, maxiter=1000, restart = 1000)
    solver_time = time.time() - solver_time
    print "The linear system was solved in {0} iterations in {:5.2f} seconds".format(it_count, solver_time)

    # The following grid function stores the computed boundary data of the total field.
    p1, p2 = np.split(x[:(dirichl_space_in+neumann_space_in)], 2) 
    total_dirichl_in = bempp.api.GridFunction(dirichl_space_in, coefficients=p1)
    total_neumann_in = bempp.api.GridFunction(neumann_space_in, coefficients=p2)

    # Calculate potentials in the coordinates of the atoms
    slp_ev = bempp.api.operators.potential.laplace.single_layer(neumann_space_in, x_q.transpose())
    dlp_ev = bempp.api.operators.potential.laplace.double_layer(dirichl_space_in, x_q.transpose())

    # Evaluate potential at charges position & total dissolution energy
    phi_q = (1./ep)*slp_ev*total_neumann_in - dlp_ev*total_dirichl_in
    total_energy = 2*np.pi*332.064*np.sum(q*phi_q).real

    total_time = time.time() - total_time
    print "Total time: {:5.2f}".format(total_time)
    print "Total dissolution Energy: {:8.3f}".format(total_energy)

    if info==True:
        info_dict = {}
        info_dict['mol_name'] = mol_name
        info_dict['mesh_density'] = mesh_density
        info_dict['n_of_elements'] = dirichl_space.global_dof_count
        info_dict['formulation'] = formulation
        info_dict['energy'] = total_energy
        info_dict['assemble_time'] = matrix_time
        info_dict['solver_time'] = solver_time
        info_dict['total_time'] = total_time
        info_dict['iterations'] = it_count
        info_dict['probe_rad'] = probe_radius

        PBL.save_log(mol_name, info_dict)

    return total_energy

mol_name = '1crn'

densities = np.array([1.])
energy = np.zeros(len(densities))

for i in range(len(densities)):
    array_it, array_frame, it_count = np.array([]), np.array([]), 0
    energy[i] = PBL_stern(mol_name, densities[i], info=True)

