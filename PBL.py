import general_functions as gf
import numpy as np
import bempp.api, inspect, time, os, sys

from scipy.sparse.linalg import gmres
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

def solvation_energy(mol_name, mesh_density, 
					 ep_in=4., ep_ex=80., kappa=0.125, 
					 solver_tol=1e-3, formulation='direct', 
					 info=False, phi_info=False ):
	'''
	mol_name  : molecule name to call .pqr & .msh files from its directories
	mesh_density : mesh density from msms/nano cration
	ep_in : interior electrostatic permittivity
	ep_ex : exterior electrostatic permittivity
	kappa : kappa constant for Poisson-Boltzmann field
	'''

	mesh_directory = 'Molecule/' + mol_name + '/mesh'
	mesh_file = '{}/{}_d{:04.1f}.msh'.format(mesh_directory, mol_name, mesh_density)

	total_time = time.time()

	# Import .msh file
	grid = bempp.api.import_grid(mesh_file)

	# Define potential and derivate spaces
	dirichl_space = bempp.api.function_space(grid, "DP", 0)
	neumann_space = bempp.api.function_space(grid, "DP", 0)

	print "\nNew Calculation for {} using {} formulation".format(mol_name, formulation)
	print "Mesh file: {}".format(mesh_file)
	print "Mesh density: {:5.2f}".format(mesh_density)
	print "Number of elements: {0}".format(dirichl_space.global_dof_count)

	q, x_q = gf.read_pqr(mol_name)

	print "Assembling Matrix..."
	matrix_time = time.time()

	if formulation=='direct':
		A, rhs = direct_formulation(dirichl_space, neumann_space, q, x_q, ep_in, ep_ex, kappa)
	elif formulation=='juffer':
		A, rhs = juffer_formulation(dirichl_space, neumann_space, q, x_q, ep_in, ep_ex, kappa)

	matrix_time = time.time() - matrix_time
	print "Assamble time: {:5.2f}".format(matrix_time)

	# print "Preconditioning"
	# precond_time = time.time()
	# from scipy.sparse import block_diag
	# A_prec = gf.inverse_block_diagonals(bempp.api.as_matrix(A).real)
	# precond_time = time.time() - precond_time

	solver_time = time.time()
	global array_it, array_frame, it_count
	array_it, array_frame, it_count = np.array([]), np.array([]), 0
	x, _ = gmres(A, rhs, callback=iteration_counter, tol=solver_tol, maxiter=500, restart = 1000)
	solver_time = time.time() - solver_time
	print "The linear system was solved in {:5.2f} seconds and {} iteration".format(solver_time, it_count)

	# Evaluate potential at charges position & total dissolution energy
	phi_q = gf.charges_potential(x, x_q, dirichl_space, neumann_space)
	total_energy = 2*np.pi*332.064*np.sum(q*phi_q)

	print "Total dissolution Energy: {:8.3f}".format(total_energy)

	total_time = time.time() - total_time
	print "Total time: {:5.2f}".format(total_time)

	if phi_info:
		if not os.path.exists('log/'):
			os.makedirs('log/')
		phi_data = open('log/phi_{}_{}_d{:04.1f}'.format(mol_name, formulation, mesh_density), 'w')
		for boundary_element in x:
			phi_data.write(str(boundary_element.real) + '\n')

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

		gf.save_log(mol_name, info_dict)

	return total_energy

def direct_formulation(dirichl_space, neumann_space, q, x_q, ep_in, ep_ex, kappa):

	# Functions to proyect the carges potential to the boundary
	def green_func(x, n, domain_index, result):
		result[:] = np.sum(q/np.linalg.norm( x - x_q, axis=1 ))/(4.*np.pi*ep_in)

	charged_grid_fun = bempp.api.GridFunction(dirichl_space, fun=green_func)
	rhs = np.concatenate([charged_grid_fun.coefficients, np.zeros(neumann_space.global_dof_count)])

	# Define Single & Double Layer Operators
	identity = sparse.identity(dirichl_space, dirichl_space, dirichl_space)
	slp_in = laplace.single_layer(neumann_space, dirichl_space, dirichl_space)
	dlp_in = laplace.double_layer(dirichl_space, dirichl_space, dirichl_space)
	slp_ex = modified_helmholtz.single_layer(neumann_space, dirichl_space, dirichl_space, kappa)
	dlp_ex = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dirichl_space, kappa)

	blocked = bempp.api.BlockedOperator(2, 2)
	blocked[0, 0] = 0.5*identity + dlp_in
	blocked[0, 1] = -slp_in
	blocked[1, 0] = 0.5*identity - dlp_ex
	blocked[1, 1] = (ep_in/ep_ex)*slp_ex
	A = blocked.strong_form()

	tree = bempp.api.hmatrix_interface.block_cluster_tree(dlp_in.weak_form)
	tree.plot()

	return A, rhs

def juffer_formulation(dirichl_space, neumann_space, q, x_q, ep_in, ep_ex, kappa):
	phi_id = sparse.identity(dirichl_space, dirichl_space, dirichl_space)
	dph_id = sparse.identity(neumann_space, neumann_space, neumann_space)
	ep = ep_ex/ep_in

	dF = laplace.double_layer(dirichl_space, dirichl_space, dirichl_space)
	dP = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dirichl_space, kappa)
	L1 = ep*dP - dF

	F = laplace.single_layer(neumann_space, dirichl_space, dirichl_space)
	P = modified_helmholtz.single_layer(neumann_space, dirichl_space, dirichl_space, kappa)
	L2 = F - P

	ddF = laplace.hypersingular(dirichl_space, neumann_space, neumann_space)
	ddP = modified_helmholtz.hypersingular(dirichl_space, neumann_space, neumann_space, kappa)
	L3 = ddP - ddF  # Cambio de signo por definicion de bempp

	dF0 = laplace.adjoint_double_layer(neumann_space, neumann_space, neumann_space)
	dP0 = modified_helmholtz.adjoint_double_layer(neumann_space, neumann_space, neumann_space, kappa)
	L4 = dF0 - (1./ep)*dP0

	blocked = bempp.api.BlockedOperator(2, 2)
	blocked[0, 0] = 0.5*(1. + ep)*phi_id - L1
	blocked[0, 1] = -L2
	blocked[1, 0] = -L3
	blocked[1, 1] = 0.5*(1. + 1./ep)*dph_id - L4
	A = blocked.strong_form()

	def d_green_func(x, n, domain_index, result):
		const = -1./(4.*np.pi*ep_in)
		result[:] = const*np.sum(q*np.dot( x - x_q, n )/(np.linalg.norm( x - x_q, axis=1 )**3))

	def green_func(x, n, domain_index, result):
		result[:] = np.sum(q/np.linalg.norm( x - x_q, axis=1 ))/(4.*np.pi*ep_in)

	rhs_1 = bempp.api.GridFunction(dirichl_space, fun=green_func)
	rhs_2 = bempp.api.GridFunction(dirichl_space, fun=d_green_func)
	rhs = np.concatenate([rhs_1.coefficients, rhs_2.coefficients])

	return A, rhs

array_it, array_frame, it_count = np.array([]), np.array([]), 0
def iteration_counter(x):
    global array_it, array_frame, it_count
    it_count += 1
    frame, array_it = inspect.currentframe().f_back, np.append(array_it, it_count)
    array_frame = np.append(array_frame, frame.f_locals["resid"])
    sys.stdout.flush()
    sys.stdout.write("Iter: {0} - Error: {1:.2E}    \r".format(it_count, frame.f_locals["resid"]))
