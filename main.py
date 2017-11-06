import numpy as np
import bempp.api
import PBL

mol_name = '5pti'
densities = np.array([ 5.7 ])

energy_d = np.zeros(len(densities))
energy_j = np.zeros(len(densities))
times_d = np.zeros((len(densities), 3))
times_j = np.zeros((len(densities), 3))
n_boundary_elements = np.zeros(len(densities))

for i in range(len(densities)):
	energy_d1 = PBL.solvation_energy(mol_name, densities[i], info=True)
   	energy_d2 = PBL.solvation_energy(mol_name, densities[i], 
   									   formulation='juffer', info=True)


directory = 'Molecule/{}/'.format(mol_name)
log_file = open(directory + mol_name + '_log', 'r').read().split('\n')
dicts = [eval(line) for line in log_file if len(line)!= 0]

n_boundary_elements = np.array([l['n_of_elements']
								 for l in dicts if l['formulation']=='direct'])

energy_d = np.array([l['energy'] for l in dicts if l['formulation']=='direct'])
energy_j = np.array([l['energy'] for l in dicts if l['formulation']=='juffer'])

times_d = np.array([[l['assemble_time'], l['solver_time'], l['total_time']]
							 for l in dicts if l['formulation']=='direct'])
times_j = np.array([[l['assemble_time'], l['solver_time'], l['total_time']]
							 for l in dicts if l['formulation']=='juffer'])

print 'N of elements: ', n_boundary_elements
print 'Direct form energy: ', energy_d
print 'Juffer form energy:', energy_j

rich_energy_d, _, _ = PBL.richardson_extrapolation(energy_d, n_boundary_elements)
rich_energy_j, _, _ = PBL.richardson_extrapolation(energy_j, n_boundary_elements)
r_solution_d = rich_energy_d*np.ones(len(n_boundary_elements))
r_solution_j = rich_energy_j*np.ones(len(n_boundary_elements))


import matplotlib.pyplot as plt
plt.switch_backend('agg')

params = {'figure.figsize': (7, 6),
		  'axes.titlesize':  12,
          'axes.labelsize':  10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
		  'legend.fontsize': 10,
		  'lines.linewidth': .5,
		  'legend.loc': 'upper left'}
plt.rcParams.update(params)

energy_rslt = plt.figure().add_subplot(111)
energy_rslt.plot(n_boundary_elements, energy_d, marker='o', label='direct', color='k')
energy_rslt.plot(n_boundary_elements, energy_j, marker='s', label='Juffer', color='k')

energy_rslt.plot(n_boundary_elements, r_solution_d, 'r--', marker='o', color='k')
energy_rslt.plot(n_boundary_elements, r_solution_j, 'r--', marker='s', color='k')

energy_rslt.set_title('Mesh Convergence')
energy_rslt.set_xlabel('N of elements')
energy_rslt.set_ylabel('G  [kcal]')
energy_rslt.legend()
plt.savefig(directory + 'Energy.png')


time_general = plt.figure().add_subplot(111)
time_general.plot(n_boundary_elements, times_d[:,0], 'r--', marker='o', label='direct assemble', color='k')
time_general.plot(n_boundary_elements, times_d[:,1], 'r:', marker='o',label='direct solver', color='k')
time_general.plot(n_boundary_elements, times_d[:,2], marker='o', label='direct total', color='k')

time_general.plot(n_boundary_elements, times_j[:,0], 'r--', marker='s',label='juffer assemble', color='k')
time_general.plot(n_boundary_elements, times_j[:,1], 'r:', marker='s', label='juffer solver', color='k')
time_general.plot(n_boundary_elements, times_j[:,2], marker='s', label='juffer total', color='k')

time_general.set_title('Time Consuming')
time_general.set_xlabel('N of elements')
time_general.set_ylabel('Time [s]')
time_general.legend()
plt.savefig(directory + 'Time.png')


# time_solving = plt.figure().add_subplot(111)
# time_solving.plot(n_boundary_elements, times_d[:,1], label='direct')
# time_solving.plot(n_boundary_elements, times_j[:,1], label='Juffer')

# time_solving.set_title('Solver Time')
# time_solving.set_xlabel('N of elements')
# time_solving.set_ylabel('Time [s]')
# time_solving.legend()
# plt.savefig(directory + 'SolverTime.png')

# time_assemble = plt.figure().add_subplot(111)
# time_assemble.plot(n_boundary_elements, times_d[:,0], label='direct')
# time_assemble.plot(n_boundary_elements, times_j[:,0], label='Juffer')

# time_assemble.set_title('Assemble Time')
# time_assemble.set_xlabel('N of elements')
# time_assemble.set_ylabel('Time [s]')
# time_assemble.legend()
# plt.savefig(directory + 'AseembleTime.png')
