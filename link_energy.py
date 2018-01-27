import PBL, PBL_stern

dens = np.array([ 2. ])#, 4., 8., 16., 32. ])

link_energy = []

for dd in dens:
	energy_p = PBL.solvation_energy('peptide2', dd, info=True)
	energy_r = PBL.solvation_energy('rna', dd, info=True)

	energy_c = PBL.solvation_energy('complex', dd, info=True)

	link_energy.append(energy_c - (energy_p + energy_r))

exit()

for dd in dens:
	energy_p = PBL_stern.solvation_energy('peptide2', dd, info=True)
	energy_r = PBL_stern.solvation_energy('rna', dd, info=True)

	energy_c = PBL_stern.solvation_energy('complex', dd, info=True)

	link_energy.append(energy_c - (energy_p + energy_r))
