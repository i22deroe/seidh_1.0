"""Seidh molecule adapt function prior to SD energy scoring function

Based on SwarmDock energy scoring function
Reference: SwarmDock and the use of Normal Models in Protein-Protein Docking
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996808/

SeidhAdapt takes exactly 1 PDB object and returns a list containing:
-The atom list from the PDB object
-The electronic charges 
-The van der waals energies
-The van der waals radii
-The covalent radii
"""

import numpy as np
from sd.data.amber import amber_types, masses, charges
import sd.data.vdw as vdw
import sd.data.cov as cov

def SeidhAdapt(molecule):
	atoms = molecule.get_atoms()
	atoms_=list() #Since the get_atoms() method provides a generator, another variable is created as a means to store the actual values.
	#Atom properties
	for atom in atoms:
		if  atom.get_parent().id[0] != ' ': #Hetero atoms are ignored
			continue
		res_name = atom.get_parent().get_resname()
		if res_name == "HOH": #Water is ignored
			continue
		if res_name == "HIS":
			res_name = 'HID'
		atom_name=atom.get_name()
		if atom_name == "OXT":
			atom_name = "O"
		atom_id = "%s-%s" % (res_name, atom_name)
		atom.amber_type = amber_types[atom_id]
		atom.charge = charges[atom_id]
		atom.mass = masses[atom.amber_type]
		atom.vdw_energy = vdw.vdw_energy[atom.amber_type]
		atom.vdw_radius = vdw.vdw_radii[atom.amber_type]
		#atom.cov_radius = cov.cov_radii[atom.amber_type]
		atoms_.append(atom)
	#Model properties
	coordinates = np.array([atom.get_vector() for atom in atoms_]) #Atom coordinates, index 0
	elec_charges = np.array([atom.charge for atom in atoms_]) #Electric charges, index 1
	vdw_energies = np.array([atom.vdw_energy for atom in atoms_]) #Van der waals energies, index 2
	vdw_radii = np.array([atom.vdw_radius for atom in atoms_]) #Van der waals radii, index 3
	#cov_radii = np.array([atom.cov_radius for atom in atoms_]) #Covalent radii, index 4
	return(coordinates,elec_charges,vdw_energies,vdw_radii) 