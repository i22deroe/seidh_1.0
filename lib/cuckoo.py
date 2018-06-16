'''
	Implementation of the cuckoo algorithm for seidh.

	The cuckoo function is called via seidh main function. It requires a fasta sequence,
	the binding point, a SeidhAdapt-ed protein and the number of procs to operate with.

	Please note that parting from Biopython's PDB-structured molecules none of these functions
	should be normally called.
	
	Glossary:
		peptide: the structure being currently predicted
		protein: the reference, target structure to which the peptide shall bind
		aminoacid: a single fasta character
		cuckoo: the current aminoacid being modified in order to bind to the peptide
		flock: a list of cuckoos
		alt_: prefix for an altered aminoacid which will be compared to the original one
		adapted_: prefix for the seidh-working structure containing the coordinates, radii, energies, etc. of a peptide
		_energy: a suffix for an aminoacid, flock or cuckoo which will be used to provide ranking criteria

		It requires the following external packages:
		-Biopython
		-PeptideBuilder
		-Nummpy
		-Multiprocessing
'''


from Bio.PDB import *
from Bio.SeqUtils import *
from sd.driver import SeidhAdapt
from sd.scoring import seidh_calculate_energy
from distance import VdW_minima,inner_stability
from levy_distribution import levy_distro as levy
from PeptideBuilder import Geometry, PeptideBuilder
from copy import deepcopy
from props.props import proteins,peptides,predictions,entities
import numpy
import multiprocessing as mp
import random

def new_cuckoo(aminoacid):
	'''
		A cuckoo is born
	'''
	cuckoo=Geometry.geometry(aminoacid)
	cuckoo.phi=random.uniform(-180.0,180.0)
	cuckoo.psi_im1=random.uniform(-180.0,180.0)
	return cuckoo
	
def modify_cuckoo(cuckoo):
	'''
		Lets a cuckoo perform a levy flight and change its angles
	'''
	phi=cuckoo.phi
	psi_im1=cuckoo.psi_im1
	cuckoo.phi=cuckoo.phi+levy(0)
	if cuckoo.phi > 180. or cuckoo.phi < -180.:	#If phi turns into a non-valid value, it is substitued by the original value.
		cuckoo.phi = phi
	cuckoo.psi_im1=cuckoo.psi_im1+levy(0)
	if cuckoo.psi_im1 > 180. or cuckoo.phi < -180: #If psi_im1 turns into a non-valid value, it is substitued by the original value.
		cuckoo_psi_im1 = psi_im1
	return cuckoo

def new_flock(peptide,aminoacid,protein):
	'''
		Returns a list of random born cuckoos from the aminoacid in line; as well as their adapted versions for seidh
	'''
	flock=list()
	adapted_flock=list()
	i=0
	while i<MAX_AA: #MAX_AA cuckoos are born!
		tmp=deepcopy(peptide)
		egg=new_cuckoo(aminoacid)
		semi_egg=PeptideBuilder.add_residue(tmp,egg)
		adapted_egg=SeidhAdapt(semi_egg)
		attempt=0
		while not VdW_minima(adapted_egg,protein) or not inner_stability(semi_egg,adapted_egg): #Collision detection
			if attempt>10: #Just 10 tries
				return False,False
			tmp=deepcopy(peptide)
			egg=new_cuckoo(aminoacid)
			semi_egg=PeptideBuilder.add_residue(tmp,egg)
			adapted_egg=SeidhAdapt(semi_egg)
			attempt+=1
		flock.append(egg) #A cuckoo egg just hatched!
		adapted_flock.append(adapted_egg)
		i+=1
	return flock,adapted_flock

def flock_migration(flock,adapted_flock,protein,peptide):
	'''
		Returns a list of cuckoos with alterations if proved to be better than the original cuckoos
	'''
	flock_energy=list()
	i=0
	for cuckoo in flock: #Cuckoos are flying!
		tmp=deepcopy(peptide)
		cuckoo_energy=seidh_calculate_energy(protein,adapted_flock[i]) #time consuming function
		alt_cuckoo=deepcopy(cuckoo)
		alt_cuckoo=modify_cuckoo(alt_cuckoo) #Now we have two cuckoos: cuckoo and its possible future self, the alt_cuckoo
		alt_peptide=PeptideBuilder.add_residue(tmp,alt_cuckoo)
		adapted_alt_peptide=SeidhAdapt(alt_peptide)
		if not VdW_minima(protein,adapted_alt_peptide) or not inner_stability(alt_peptide,adapted_alt_peptide): #The altered cuckoo collides, and hence is discarded
			flock_energy.append(cuckoo_energy)
		else: #altered cuckoo is ok!
			alt_cuckoo_energy=seidh_calculate_energy(protein,adapted_alt_peptide) #time consuming function
			if alt_cuckoo_energy>cuckoo_energy:
				cuckoo=alt_cuckoo #alt_cuckoo becomes a real cuckoo!
				flock_energy.append(alt_cuckoo_energy)
			else:
				flock_energy.append(cuckoo_energy)
		i+=1
	return flock_energy,flock
	
def cuckoo_migration(cuckoo,adapted_cuckoo,protein,peptide):
	'''
		Returns a single cuckoo with alterations, if proved to be better than the original cuckoo
	'''
	tmp=deepcopy(peptide)
	
	cuckoo_energy=seidh_calculate_energy(protein,adapted_cuckoo) #time consuming function
	
	alt_cuckoo=deepcopy(cuckoo)
	alt_cuckoo=modify_cuckoo(alt_cuckoo) #Now we have two cuckoos: cuckoo and its possible future self, the alt_cuckoo
	alt_peptide=(PeptideBuilder.add_residue(tmp,alt_cuckoo))
	adapted_alt_peptide=SeidhAdapt(alt_peptide)
	
	if not VdW_minima(protein,adapted_alt_peptide) or not inner_stability(alt_peptide,adapted_alt_peptide): #In case the altered cuckoo is too close, avoids a time consuming function
		return cuckoo_energy,cuckoo
	else: #Altered cuckoo is ok!
		alt_cuckoo_energy=seidh_calculate_energy(protein,adapted_alt_peptide) #time consuming function
		if alt_cuckoo_energy>cuckoo_energy:
			return alt_cuckoo_energy,alt_cuckoo
		else:
			return cuckoo_energy,cuckoo
			
def multi_migration(flock,adapted_flock,protein,peptide, NUM_PROC):
	pool = mp.Pool(processes=NUM_PROC)
	results = [pool.apply_async(cuckoo_migration, args=(cuckoo,adapted_cuckoo,protein,peptide)) for cuckoo,adapted_cuckoo in zip(flock,adapted_flock)]
	results = [p.get() for p in results]
	cuckoo_list,energy_list = zip(*results)
	'''
	cuckoo_list=list()
	energy_list=list()
	
	for ele in results:
		cuckoo_list.append(ele[0])
		energy_list.append(ele[1])
	'''
	pool.close()
	
	return cuckoo_list,energy_list
	
def rank(flock_energy,flock,diff=1):
	ranked=[cuckoo for _,cuckoo in sorted(zip(flock_energy,flock),reverse=True)]
	if diff==0:
		return ranked[:-10]
	else: 
		return ranked

def cuckoo (aminoacid_sequence, peptide, protein, NUM_PROC):
	flock=list()
	adapted_flock=list()
	i=0
	for aminoacid in aminoacid_sequence:
		print "%.2f%% complete" % (i/float(len(aminoacid_sequence))*100)
		flock,adapted_flock=new_flock(peptide,aminoacid,protein) #First generation of cuckoos are born!
		if flock is None or flock==0:
			return False
		while len(flock)>10:
			if NUM_PROC == 1:
				flock_energy,flock=flock_migration(flock,adapted_flock,protein,peptide)
			else:
				flock_energy,flock=multi_migration(flock,adapted_flock,protein,peptide,NUM_PROC)
			flock=rank(flock_energy,flock,0)
			adapted_flock=rank(flock_energy,adapted_flock,0)
		peptide=PeptideBuilder.add_residue(peptide,flock.pop(0))
		i+=1
	return peptide
