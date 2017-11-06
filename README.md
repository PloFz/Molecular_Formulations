# Molecular-Calculations
Code to calculate the electrostatic potential and dissolution energy in a molecule immerse in a dielectric medium. Solved by a Boundary Element Method and implemented in Python and BEM++.

## Folders:
### src_files
	*.pqr : Position, Charge & Raduis data (obtained by pdb2pqr)
	*.pbd : Original data file from Protein Data Bank
### Mesh_(mol name)
	Mesh files in format: mol-name_dentisty_(stern) (obtained by msms)
	*.vert & *.face : msms output, stores the vertex and faces for each mesh
	*.msh : GMESH file, obtained from .vert & .face via GridFactory() function
### log
	Record basic calculations data, file name in foramt: mol-name_density_date:hour
## Programs
### PLB.py 
	Function to calculate Dissolution Energy
		input : .pqr and .msh
		output : Returns the Energy in [kcal/mol] and print a log file

## PLB_stern.py
	Same PLB with an additional Stern Layer
