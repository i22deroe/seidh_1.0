""" 
	This file gathers a few functions that are used by Seidh, or that will be used in
	future updates.
	
	It requires the following external packages:
		-Biopython
		-PeptideBuilder
	
"""

#Libs
import os
from Bio.PDB import *
from Bio.SeqUtils import *
from PeptideBuilder import Geometry, PeptideBuilder
from props.props import proteins,peptides,predictions,entities
import random

#Constants
fasta_ambiguous=('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','X','Z')
fasta_unambiguous=('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y')
ambiguity_dictionary={'B':('D','N'),'J':('L','I'),'Z':('E','Q'),'X':fasta_ambiguous}

#Functions
def is_fasta(fasta):
	'''
		This function recognizes whether a fasta sequence is composed of fasta symbols or not.
	'''
	for symbol in fasta:
		if symbol not in fasta_ambiguous:
			return False
	return True

def fasta_ambiguity(fasta):
	'''
		This function returns an unambiguous fasta sequence. Since the fasta code may contain wildcards, it is precise to determine
			an unambiguous symbol prior to processing.
	'''
	if is_fasta(fasta):
		new_seq=list()
		for symbol in fasta:
			if symbol in ambiguity_dictionary.keys():
				rnd_val=random.randint(0,len(fasta_unambiguous))
				new_seq.append(ambiguity_dictionary[symbol][rnd_val%len(ambiguity_dictionary[symbol])])
			else:
				new_seq.append(symbol)
		return new_seq
	else:
		return False
		
def open_pdb_protein(file):
	'''
		Open a PDB file containing a protein. 
	'''
	if not os.path.isfile(proteins+file):
		return False
	else:
		parser=PDBParser()
		return parser.get_structure(file,proteins+file)

def open_pdb_entity(file,chain):
	'''
		Open a pdb entity provided a file name and a chain name, returning the structure for the sole chain.
	'''
	if not os.path.isfile(entities+"pdb"+file+".ent"):
		print entities+"pdb"+file+".ent"
		return False
	else:
		parser=PDBParser()
		structure=parser.get_structure(file,entities+"pdb"+file)
		return structure[0][chain]

def open_pdb_peptide(file):
	'''
		Open a PDB file containing an amino acid or peptide. 
	'''
	if not os.path.isfile(peptides+file):
		return False
	else:
		parser=PDBParser()
		structure=parser.get_structure(file,peptides+file)
		return structure

def open_fasta_peptide(file):
	'''
		Open the fasta sequence of a peptide if it's within the peptides directory or try instead to read as provided.
	'''
	if not os.path.isfile(peptides+file):
		return fasta_ambiguity(file)
	else:
		fasta_file = open(peptides+file,'r')
		sequence = fasta_file.read()
		fasta_file.close()
		return fasta_ambiguity(sequence)

def cut_peptide(peptide,fasta):
	'''
		Returns the remaining fasta sequence to append to the given peptide.
	'''
	if is_fasta(fasta):
		sequence=list()
		for residue in peptide.get_residues():
			sequence.append(residue.get_resname())
		sequence=seq1(sequence)
		return sequence[len(fasta):]
	else:
		return False
		
def save_structure(structure,name):
	'''
		Save a PDB structure to a PDB file
	'''
	io=PDBIO()
	io.set_structure(structure[0]['A'])
	io.save(predictions+name+'.pdb')
	return True
