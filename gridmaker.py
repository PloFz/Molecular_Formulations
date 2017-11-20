import numpy as np, bempp.api, os

def pqr2mesh(mol_name, density=3.,
			 probe_radius=1.4, min_area=1e-5,
			 stern=False, stern_radius=1., program='nano'):
	'''
	This function returns a bempp grid object and creates a .msh files in Mesh/ directory
	from a .pqr file, the following inputs can me modified:
		mol_name 	 : molecule name
		pqr_file 	 : name and directory of the .pqr file
		probe_radius : probe atom radius 
		density 	 : grid density 
		min_area 	 : minimum boundary elements size to avoid singular matrix
		stern 		 : stern layer radius
	'''

	# Directories
	mol_directory = 'Molecule/' + mol_name
	mesh_directory = mol_directory + '/mesh'
	geom_directory = mol_directory + '/geometry'
	nano_directory = 'ExternalPrograms/NanoShaper/'
	temp_directory = 'mesh_temp/'

	if not os.path.exists(temp_directory):
		os.makedirs(temp_directory)
	if not os.path.exists(mesh_directory):
		os.makedirs(mesh_directory)
	if not os.path.exists(geom_directory):
		os.makedirs(geom_directory)

	# Files
	pqr_file = mol_directory + '/' + mol_name + '.pqr'
	xyzr_file = temp_directory + '/' + mol_name + '.xyzr'
	mesh_name = '{}/{}_d{:04.1f}'.format(mesh_directory, mol_name, density)
	geom_name = '{}/{}_d{:04.1f}'.format(geom_directory, mol_name, density)

	# Create .xyzr
	atoms_file = open(pqr_file, 'r').read().split('\n')
	xyrz_data = open(xyzr_file, 'w')
	for line in atoms_file:
			line = line.split()
			if len(line)==0 or line[0]!='ATOM':	continue		# do not read empty lines
			if stern:
				line[9] = str(float(line[9]) + stern_radius)	# add 'stern' to radius
			xyrz_data.write(line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[9]+"\n")
	xyrz_data.close()

	if stern:
		mesh_name += "_strn-pr{:04.1f}".format(stern_radius)
		geom_name += "_strn-pr{:04.1f}".format(stern_radius)

	if program == 'msms':
		# Write msms command to create .vert & .face files
		msms, mode = "~/.msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1 ", "-no_header "
		prob_rad, dens_msh = " -probe_radius " + str(probe_radius), " -density " + str(density)
		os.system( msms + mode + "-if " + xyzr_file + " -of "+ geom_name)

	# Execute NanoShaper
	if program == 'nano':
		#os.system('mv ' + xyzr_file + ' ExternalPrograms/NanoShaper/')
		config_file = open('ExternalPrograms/NanoShaper/config', 'r')
		new_config = open( temp_directory + 'surfaceConfiguration.prm', 'w')

		for line in config_file:
			if 'XYZR_FileName' in line:
				line = 'XYZR_FileNae = ' + mol_name + '.xyzr \n'
			elif 'Grid_scale' in line:
				line = 'Grid_scale = {:04.1f} \n'.format(density)
			elif 'Probe_Radius' in line:
				line = 'Probe_Radius = {:03.1f} \n'.format(probe_radius)
			new_config.write(line)
		
		new_config.close()
		config_file.close()

		os.chdir(temp_directory)
		os.system('ls')
		os.system('head surfaceConfiguration.prm')
		os.system('./../ExternalPrograms/NanoShaper/NanoShaper')
		os.chdir('..')

		os.system('mv ' + temp_directory + '*.vert ' + geom_name + '.vert')
		os.system('mv ' + temp_directory + '*.face ' + geom_name + '.face')
		os.system('rm -r ' + temp_directory)

	vertx_file = open(geom_name + ".vert", 'r').read().split('\n')
	faces_file = open(geom_name + ".face", 'r').read().split('\n')

	# Counters for small triangles, and total ignored area
	xcount, atotal, a_excl = 0, 0., 0.
	vertex = np.empty((0,3))

	# Create the grid with the factory method
	factory = bempp.api.grid.GridFactory()
	for line in vertx_file:
	    line = line.split()
	    if len(line) != 9: continue
	    vertex = np.vstack(( vertex, np.array([line[0:3]]).astype(float) ))
	    factory.insert_vertex(vertex[-1])

	# Grid assamble, exclude elements < min_area
	for line in faces_file:
		line = line.split()
		if len(line) != 5: continue
		A, B, C, _, _ = np.array(line).astype(int)
		side1, side2  = vertex[B-1] - vertex[A-1], vertex[C-1] - vertex[A-1]
		face_area = 0.5*np.linalg.norm(np.cross(side1, side2))
		atotal += face_area
		if face_area > min_area:
			factory.insert_element([A-1, B-1, C-1])
		else:
			xcount += 1			 # excluded element counter not used
			a_excl += face_area  # total area excluded not used

	grid = factory.finalize()
	
	bempp.api.export(grid=grid, file_name=mesh_name + '.msh')
	return grid

# Create a .msh file from a .pqr file
mol_name = "5pti"
mol_directory = "Molecule/" + mol_name + "/"
pqr_file_name = mol_directory + mol_name + ".pqr"

# If its a new molecule, create a directory and a .pqr file from a .pdb with pdb2pqr
if not os.path.exists(mol_directory):
	os.makedirs(mol_directory)
if not os.path.exists(pqr_file_name):
	pdb2pqr, method = "~/.pdb2pqr-linux-bin64-2.1.0/pdb2pqr ", "--ff=amber "
	pdb_file_name = "PDB_files/" + mol_name + ".pdb "
	os.system( pdb2pqr + method + pdb_file_name + pqr_file_name )


#dens = [ .8, 1., 2., 2.8, 4., 5.7 ]
dens = [ 1. ]
for dd in dens:
	grid_in = pqr2mesh(mol_name, density=dd, program='msms')
	grid_in = pqr2mesh(mol_name, density=dd, stern=True, stern_radius=1.4, program='msms')
# r_st = [ 1.4 ]
# #r_st = [ .1, .2, .4, .6, .8, 1.2, 1.5, 2., 3., 4. ]
# for rr in r_st:
# 	grid_ex = pqr2mesh(mol_name, density=5.7, stern=True, stern_radius=rr, program='msms')

#grid_in.plot()

print 'Finished'
