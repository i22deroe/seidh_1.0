"""Seidh's energy calculation function

Based on SwarmDock energy scoring function
Reference: SwarmDock and the use of Normal Models in Protein-Protein Docking
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996808/

Seidh's energy calculation function takes exactly 2 Seidh adapted PDB objects and returns a float.

It can be called the following way: seidh_calculate_energy(SeidhAdapt(protein),SeidhAdapt(peptide))

"""


from numpy import *

EPSILON=4.0
FACTOR=332.0
CUTON=7.0
CUTOFF=9.0
CUTON2=CUTON*CUTON
CUTOFF2=CUTOFF*CUTOFF
VDW_CUTOFF=5000000.0

def seidh_calculate_energy(rec,lig):

	energy=0.
	
	if rec is not None and lig is not None:
		rec_len=len(rec[0])
		lig_len=len(lig[0])

		total_elec=0.
		atom_elec=0.
		total_vdw=0.
		
		i=0
		while i<rec_len:
			atom_vdw=0.
			j=0
			while j<lig_len:
				x=rec[0][i][0] - lig[0][j][0]
				y=rec[0][i][1] - lig[0][j][1]
				z=rec[0][i][2] - lig[0][j][2]
				distance= x*x+y*y+z*z
				if distance<CUTOFF2:
					atom_elec=rec[1][i]*lig[1][j]/distance
					atom_elec*=FACTOR/EPSILON
					vdw_energy=sqrt(rec[2][i]*lig[2][j])
					vdw_radius=rec[3][i]*lig[3][j]
					p6=pow(vdw_radius,6)/pow(distance,3)
					k=vdw_energy*(p6*p6-2.0*p6)
					atom_vdw+=k
					if atom_vdw > VDW_CUTOFF:
						atom_vdw=VDW_CUTOFF
					if distance < CUTON2:
						energy+=atom_elec+atom_vdw
					else:
						energy+=(atom_elec+atom_vdw)*( (CUTOFF2 - distance)*(CUTOFF2 - distance) * (CUTOFF2 + 2.*distance - 3.0*CUTON2) / ((CUTOFF2-CUTON2)*(CUTOFF2-CUTON2)*(CUTOFF2-CUTON2)) )
				j+=1
			i+=1
	return energy*-1.
	

