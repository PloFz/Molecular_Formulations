import general_functions as gf
import gc

# dens = np.array([ .8, 1., 2., 2.8, 4., 5.7 ])

# for dd in dens:
# 	PBL.solvation_energy(mol_name, dd, info=True)
# 	gc.collect()
# 	PBL.solvation_energy(mol_name, dd, formulation='juffer', info=True)
# 	gc.collect()
# 	PBL_stern.solvation_energy(mol_name, dd, info=True)
# 	gc.collect()
# 	gf.run_pygbe('5pti', dd, 1.4, info=True)
# 	gc.collect()

# r_st = np.array([ .1, .2, .4, .6, .8, 1.2, 1.5, 2., 3., 4. ])
# for rr in r_st:
# 	PBL_stern.solvation_energy(mol_name, 4., stern_radius=rr, info=True)
# 	gc.collect()
# 	gf.run_pygbe('5pti', 4., rr, info=True)
# 	gc.collect()