import general_functions as gf, numpy as np, gc, PBL, PBL_stern

mol_name = '5pti'
# PBL_stern.solvation_energy(mol_name, 2., info=True)

dens = np.array([ 8., 16., 32. ])

for dd in dens:
	PBL.solvation_energy(mol_name, dd, info=True)
 	gc.collect()
 	PBL.solvation_energy(mol_name, dd, formulation='juffer', info=True)
 	gc.collect()
 	PBL_stern.solvation_energy(mol_name, dd, info=True)
 	gc.collect()
 	PBL_stern.solvation_energy(mol_name, dd, formulation='asc', info=True)
 	gc.collect()
	gf.run_pygbe(mol_name, dd, 1.4, info=True)
 	gc.collect()


r_st = np.array([ .4, .6, .8, 1.2, 1.5, 2., 3., 4. ])

for rr in r_st:
 	PBL_stern.solvation_energy(mol_name, 32., stern_radius=rr, info=True)
 	gc.collect()
 	PBL_stern.solvation_energy(mol_name, 32., stern_radius=rr, formulation='asc', info=True)
 	gc.collect()
 	gf.run_pygbe(mol_name, 32., rr, info=True)
 	gc.collect()
