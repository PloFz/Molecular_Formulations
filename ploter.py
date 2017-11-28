import general_functions as gf, numpy as np, bempp.api, PBL_stern, PBL, gc

mol_name = '5pti'
directory = 'Molecule/{}/'.format(mol_name)

log_file = open(directory + mol_name + '_log', 'r').read().split('\n')
dicts = [eval(line) for line in log_file if len(line)!= 0]

n_boundary_elements = [l['n_of_elements']
								 for l in dicts if l['formulation']=='direct']

times_d = np.array([[l['assemble_time'], l['solver_time'], l['total_time']]
								 for l in dicts if l['formulation']=='direct'])
energy_d = [l['energy'] for l in dicts if l['formulation']=='direct']

times_j = np.array([[l['assemble_time'], l['solver_time'], l['total_time']]
								 for l in dicts if l['formulation']=='juffer'])
energy_j = [l['energy'] for l in dicts if l['formulation']=='juffer']


n_stern_1 = [l['n_of_elements']
						 for l in dicts if l['formulation']=='stern_d' and l['stern_radius']==1.4]
energy_s1 = [l['energy'] for l in dicts if l['formulation']=='stern_d' and l['stern_radius']==1.4]
times_s1 = np.array([[l['assemble_time'], l['solver_time'], l['total_time']]
						 for l in dicts if l['formulation']=='stern_d' and l['stern_radius']==1.4])

n_pygbe_1 = [l['n_of_elements']
						 for l in dicts if l['formulation']=='PyGBe' and l['stern_radius']==1.4]
energy_pg = [l['energy'] for l in dicts if l['formulation']=='PyGBe' and l['stern_radius']==1.4]
times_pg = np.array([[l['solver_time'], l['total_time']]
						 for l in dicts if l['formulation']=='PyGBe' and l['stern_radius']==1.4])

print n_pygbe_1
print energy_pg

n_stern_2 = [l['n_of_elements']
						 for l in dicts if l['formulation']=='stern_d' and l['stern_radius']!=1.4]
energy_s2 = [l['energy'] for l in dicts if l['formulation']=='stern_d' and l['stern_radius']!=1.4]
radius_s2 = [l['stern_radius']
						 for l in dicts if l['formulation']=='stern_d' and l['stern_radius']!=1.4]

energy_pr = [l['energy'] for l in dicts if l['formulation']=='PyGBe' and l['stern_radius']!=1.4]
radius_pr = [l['stern_radius']
						 for l in dicts if l['formulation']=='PyGBe' and l['stern_radius']!=1.4]


rich_energy_d, _, _ = gf.richardson_extrapolation(energy_d, n_boundary_elements)
rich_energy_j, _, _ = gf.richardson_extrapolation(energy_j, n_boundary_elements)
rich_energy_s, _, _ = gf.richardson_extrapolation(energy_s1, n_stern_1)
#rich_energy_p, _, _ = gf.richardson_extrapolation(energy_pg, n_pygbe_1)

r_solution_d = rich_energy_d*np.ones(len(n_boundary_elements))
r_solution_j = rich_energy_j*np.ones(len(n_boundary_elements))
r_solution_s = rich_energy_s*np.ones(len(n_boundary_elements))
#r_solution_p = rich_energy_p*np.ones(len(n_boundary_elements))

import matplotlib.pyplot as plt
plt.switch_backend('agg')

params = {'figure.figsize':  (7, 6),
		  'axes.titlesize':  12,
          'axes.labelsize':  10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
		  'legend.fontsize': 10,
		  'lines.linewidth': .5,
		  'lines.color': 'k',
		  'legend.loc': 'upper left'}
plt.rcParams.update(params)

energy_rslt = plt.figure().add_subplot(111)
energy_rslt.plot(n_boundary_elements, energy_d, marker='o', label='direct', color='k')
energy_rslt.plot(n_boundary_elements, energy_j, marker='s', label='Juffer', color='k')
#energy_rslt.plot(n_stern_1, energy_s1, marker='^', label='Stern', color='k')
energy_rslt.plot(n_pygbe_1, energy_pg, marker='x', label='PyGBe', color='k')

energy_rslt.plot(n_boundary_elements, r_solution_d, 'r--', marker='o', color='k')
energy_rslt.plot(n_boundary_elements, r_solution_j, 'r--', marker='s', color='k')
#energy_rslt.plot(n_boundary_elements, r_solution_s, 'r--', marker='s', color='k')
#energy_rslt.plot(n_boundary_elements, r_solution_p, 'r--', marker='x', color='k')


energy_rslt.set_title('Mesh Convergence')
energy_rslt.set_xlabel('N of elements')
energy_rslt.set_ylabel(r'$\Delta$G  [kcal]')
energy_rslt.legend()
plt.savefig(directory + 'Energy.png')


time_general = plt.figure().add_subplot(111)
time_general.loglog(n_boundary_elements, times_d[:,0], 'r--', marker='o', label='direct assemble', color='k')
time_general.loglog(n_boundary_elements, times_d[:,1], 'r:', marker='o', label='direct solver', color='k')
time_general.loglog(n_boundary_elements, times_d[:,2], marker='o', label='direct total', color='k')

time_general.loglog(n_boundary_elements, times_j[:,0], 'r--', marker='s', label='juffer assemble', color='k')
time_general.loglog(n_boundary_elements, times_j[:,1], 'r:', marker='s', label='juffer solver', color='k')
time_general.loglog(n_boundary_elements, times_j[:,2], marker='s', label='juffer total', color='k')

time_general.loglog(n_stern_1, times_s1[:,0], 'r--', marker='^', label='stern assemble', color='k')
time_general.loglog(n_stern_1, times_s1[:,1], 'r:', marker='^', label='stern solver', color='k')
time_general.loglog(n_stern_1, times_s1[:,2], marker='^', label='stern total', color='k')

time_general.set_title('Time Consuming')
time_general.set_xlabel('N of elements')
time_general.set_ylabel('Time [s]')
time_general.legend(prop={'size': 13})
plt.savefig(directory + 'Time.png')


stern_fig = plt.figure()
stern_rad = stern_fig.add_subplot(111)
stern_rad.plot(radius_s2, energy_s2, marker='o', color='k')
stern_rad.plot(radius_pr, energy_pr, marker='s', color='k')

stern_rad.set_xlabel('Stern Radius')
stern_rad.set_ylabel(r'$\Delta$G  [kcal]')

plt.savefig(directory + 'Stern.png')
