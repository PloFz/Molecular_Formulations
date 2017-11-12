import numpy as np
import bempp.api
import inspect
import time
import PBL
import os

def read_grid_value(mol_name, mesh_density, formulation='direct'):
	mesh_directory = 'Molecule/' + mol_name + '/mesh'
	mesh_file = '{}/{}_d{:04.1f}.msh'.format(mesh_directory, mol_name, mesh_density)
	grid = bempp.api.import_grid(mesh_file)
	#grid = bempp.api.shapes.regular_sphere(3)

	return grid

	phi_file = 'log/phi_{}_{}_d{:04.1f}'.format(mol_name, formulation, mesh_density )
	phi_data = open(phi_file, 'r').read().strip().split('\n')
	phi_data = np.array(phi_data).astype(float)
	space = bempp.api.function_space(grid, "DP", 0)
	grid_fun = bempp.api.GridFunction(space, coefficients=phi_data)
	#return grid_fun

def read_pqr(mol_name):
	# Read charges and coordinates from .pqr
	mol_directory = 'Molecule/' + mol_name
	pqr_file = mol_directory + '/' + mol_name + '.pqr'
	pqr_file = open(pqr_file, 'r').read().split('\n')
	q, x_q = np.array([]), np.empty((0,3))
	for line in pqr_file:
		line = line.split()
		if len(line)==0 or line[0]!='ATOM': continue
		q = np.append(q, float(line[8]))
		x_q = np.vstack(( x_q, np.array(line[5:8]).astype(float) ))

	return q, x_q

def richardson_extrapolation(values, n_of_elements):
	r = n_of_elements[-1]/n_of_elements[-2]
	p = np.log((values[-3] - values[-2])/(values[-2] - values[-1]))/np.log(r)
	f_rich = values[-1] + (values[-1] - values[-2])/(r**p - 1)	#f1 fine, f2 coarse

	return f_rich, r, p

def write_dict(mol_name, _info, log_file):
	log_file.write("{{'mesh_density': {}".format(str(_info['mesh_density'])))
	log_file.write(", 'n_of_elements': {}".format(str(_info['n_of_elements'])))
	log_file.write(", 'formulation': '{}'".format(str(_info['formulation'])))
	log_file.write(", 'energy': {}".format(str(_info['energy'])))

	del _info['mesh_density']
	del _info['n_of_elements']
	del _info['formulation']
	del _info['energy']

	for key in _info:
		log_file.write(", '{}': {}".format(str(key), str(_info[key])))
	log_file.write("} \n")

def save_log(mol_name, info):
	log_path = 'Molecule/{}/{}_log'.format(mol_name, mol_name)
	
	dicts = []

	if not os.path.exists(log_path):
		log_file = open(log_path, 'w')
	else:
		log_file = open(log_path, 'r').read().split('\n')

		for line in log_file:
			if len(line)!=0:
				d = eval(line)
				if d['mesh_density']!=info['mesh_density'] or d['formulation']!=info['formulation']:
					dicts.append(d)
				elif d['formulation']=='stern_d':
					if info['formulation']=='stern_d' and d['stern_radius']!=info['stern_radius']:
						dicts.append(d)

	info_temp = info.copy()
	del info_temp['mol_name']
	dicts.append(info_temp)

	dicts = sorted(dicts, key=lambda k: k['mesh_density']) 
	dicts = sorted(dicts, key=lambda k: k['formulation']) 

	log_file = open( log_path, 'w')
	for d in dicts:
		write_dict(mol_name, d, log_file)
	log_file.close()

	
def charges_potential(x, x_q, dirichl_space, neumann_space):
	# # return evaluation: phi_k = V[d_phi] - K[phi]
	p1, p2 = np.split(x, 2)
	phi_surface  = bempp.api.GridFunction(dirichl_space, coefficients=p1)
	dphi_surface = bempp.api.GridFunction(neumann_space, coefficients=p2)

	slp_ev = bempp.api.operators.potential.laplace.single_layer(neumann_space, x_q.transpose())
	dlp_ev = bempp.api.operators.potential.laplace.double_layer(dirichl_space, x_q.transpose())
	
	phi_q = slp_ev*dphi_surface - dlp_ev*phi_surface

	return phi_q.real
