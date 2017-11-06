import numpy as np
import bempp.api
import time
import os

def PBL_stern(mol_name, mesh_density, ep_in=4., ep_ex=80., kappa=0.125):
    '''
    mol_name  : molecule name to call .pqr & .msh files from its directories
    mesh_density : mesh density from msms cration
    ep_in : interior electrostatic permittivity
    ep_ex : exterior electrostatic permittivity
    kappa : kappa constant for Poisson-Boltzmann field
    '''
    total_time = time.time()

    grid_in = bempp.api.import_grid('Mesh_{}/{}_d{:04.1f}.msh'.format(mol_name, mol_name, mesh_density))
    grid_ex = bempp.api.import_grid('Mesh_{}/{}_d{:04.1f}_stern.msh'.format(mol_name, mol_name, mesh_density))

    #grid_in.plot()
    #grid_ex.plot()

    # Define spaces for both boundary phi & dphi
    dirichl_space_in = bempp.api.function_space(grid_in, "DP", 0)
    neumann_space_in = bempp.api.function_space(grid_in, "DP", 0)
    dirichl_space_ex = bempp.api.function_space(grid_ex, "DP", 0)
    neumann_space_ex = bempp.api.function_space(grid_ex, "DP", 0)

    # number of elemnts in both grids
    n_elmts_in = dirichl_space_in.global_dof_count
    n_elmts_ex = dirichl_space_ex.global_dof_count

    # Read charges and coordinates
    pqr_file = open('src_files/' + mol_name + '.pqr', 'r').read().split('\n')
    q, x_q = np.array([]), np.empty((0,3))
    for line in pqr_file:
        line = line.split()
        if len(line)==0 or line[0]!='ATOM': continue
        q = np.append( q, float(line[8]) )
        x_q = np.vstack(( x_q, np.array(line[5:8]).astype(float) ))

    def charges_fun(x, n, domain_index, result):
        result[:] = np.sum(q/np.linalg.norm( x - x_q, axis=1 ))/(4.*np.pi*ep_in)

    # def charges_dfun(x, n, domain_index, result):
    #     const = -1./(4.*np.pi*ep_in)
    #     result[:] = const*np.sum(q*np.dot( x - x_q, n )/(np.linalg.norm( x - x_q, axis=1 )**3))

    print "Proyecting charges over surface..."
    charged_grid_fun  = bempp.api.GridFunction(dirichl_space_in, fun=charges_fun)
    # charged_grid_dfun = bempp.api.GridFunction(neumann_space_in, fun=charges_dfun)

    # charged_grid_fun.plot()
    # charged_grid_dfun.plot()

    rhs = np.concatenate([charged_grid_fun.coefficients,
			  NP.zeros(n_elmts_in),
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

    Original formulation
    blocked[1, 0] = 0.5*ident_in - dlp_in   # dlp_in = dlp_1T1
    blocked[1, 1] =  slp_in   # slp_in = slp_1T1
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
    from scipy.sparse.linalg import gmres
    x, info = gmres(A, rhs, callback=iteration_counter, tol=1e-3, maxiter=1000, restart = 1000)
    solver_time = time.time() - solver_time
    #print "The linear system was solved in {0} iterations in {:5.2f} seconds".format(it_count, solver_time)

    # The following grid function stores the computed boundary data of the total field.
    total_dirichl_in = bempp.api.GridFunction(dirichl_space_in,
    					      coefficients=x[:n_elmts_in])
    total_neumann_in = bempp.api.GridFunction(neumann_space_in,
    					      coefficients=x[n_elmts_in:2*n_elmts_in])

    # Calculate potentials in the coordinates of the atoms
    slp_ev = bempp.api.operators.potential.laplace.single_layer(neumann_space_in, x_q.transpose())
    dlp_ev = bempp.api.operators.potential.laplace.double_layer(dirichl_space_in, x_q.transpose())

    # Evaluate potential at charges position & total dissolution energy
    phi_q = slp_ev*total_neumann_in - dlp_ev*total_dirichl_in
    total_energy = 2*np.pi*332.064*np.sum(q*phi_q).real

    total_time = time.time() - total_time
    print "Total time: {:5.2f}".format(total_time)
    print "Total dissolution Energy: {:8.3f}".format(total_energy)

    # Write calculations log
    log_file = open( 'log/{}_{}_stern_d{:04.1f}'.format(time.strftime("%d%m:%H%M%S"),  mol_name, mesh_density), 'w')
    log_file.write('Molecule Name: ' + mol_name + '\n')
    log_file.write('Mesh Density: ' + str(mesh_density) + '\n')
    log_file.write('Number of elements in: ' + str(n_elmts_in) + '\n')
    log_file.write('Number of elements ex: ' + str(n_elmts_ex) + '\n')
    log_file.write('Iteration to solve: ' + str(it_count) + '\n')
    log_file.write('Assemble Time: ' + str(matrix_time) + '\n')
    log_file.write('Solver Time: ' + str(solver_time) + '\n')
    log_file.write('Total Time: ' + str(total_time) + '\n')
    log_file.write("Total dissolution Energy: {0} kcal/mol".format(total_energy))

    total_dirichl_ex = bempp.api.GridFunction(dirichl_space_ex, 
                                 coefficients=x[2*n_elmts_in:2*n_elmts_in+n_elmts_ex])
    total_neumann_ex = bempp.api.GridFunction(neumann_space_in, 
                                 coefficients=x[2*n_elmts_in+n_elmts_ex:])
    #total_dirichl_in.plot()
    #total_neumann_ex.plot()

    return total_energy

import inspect
def iteration_counter(x):
    global array_it, array_frame, it_count
    it_count += 1
    frame, array_it = inspect.currentframe().f_back, np.append(array_it, it_count)
    array_frame = np.append(array_frame, frame.f_locals["resid"])
    print "Iter: {0} - Error: {1:.2E}".format(it_count, frame.f_locals["resid"])

mol_name = '1crn'

densities = np.array([1.])
energy = np.zeros(len(densities))

for i in range(len(densities)):
    array_it, array_frame, it_count = np.array([]), np.array([]), 0
    energy[i] = PBL_stern(mol_name, densities[i])
