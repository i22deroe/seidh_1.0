import Bio.PDB
import numpy as np
from numpy.core.umath_tests import inner1d

def VdW_minima(molecule_1,molecule_2):
	'''
		This function takes two SeidhAdapt-ed molecules in order to calculate whether a collision would take place or not, depending on their VdW Radii.
	'''
	
	num_atoms_1=len(molecule_1[0])
	num_atoms_2=len(molecule_2[0])
	for i in range (num_atoms_1):
		for j in range(num_atoms_2):
			distance=np.sqrt(inner1d(list(molecule_1[0][i] - molecule_2[0][j]),list(molecule_1[0][i] - molecule_2[0][j])))
			vdw=molecule_1[3][i] + molecule_2[3][j]
			if distance < 0.3*vdw:
				return False
			j+=1
		i+=1
	return True

def molecule_contacts(peptide,protein):
	num_atoms_1=len(peptide[0])
	num_atoms_2=len(protein[0])
	contact_score=0
	
	for i in range(num_atoms_1):
		for j in range(num_atoms_2):
			distance=np.sqrt(inner1d(list(peptide[0][i] - protein[0][j]),list(peptide[0][i] - protein[0][j])))
			vdw_distance=peptide[3][i] + protein[3][j]
			if distance <= vdw_distance+0.5:
				contact_score+=1
				break
	return contact_score/float(num_atoms_1)

def inner_stability(peptide,adapted_peptide):

	distances=list()
	vdw=list()
	inventory=list()

	for residue in peptide.get_residues():
		inventory.append(0)
		for atom in residue.get_atom():
			inventory[-1]+=1

	i=0
	_i=0
	current_res=0

	while i<sum(inventory[0:len(inventory)]):
		j=sum(inventory[0:current_res+1])+1 #Avoids the first atom of the next peptide (bound atom)
		while j<len(adapted_peptide[0]):
			distances.append(np.sqrt(inner1d(list(adapted_peptide[0][i] - adapted_peptide[0][j]),list(adapted_peptide[0][i] - adapted_peptide[0][j]))))
			vdw.append(adapted_peptide[3][i] + adapted_peptide[3][j])
			j+=1
		if i==sum(inventory[0:current_res]):
			current_res+=1
		i+=1
	for element in vdw:
			if np.any(distances <= 0.4*element):
				return False
	return True

