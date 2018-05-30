import Bio.PDB
import sys
import os


def get_ca_atoms(ids_list,folder):
    ca_atoms = {}
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    for struct_id in ids_list:
        pdb_file = "%d.pdb" % struct_id
        print "Reading CA from %s" % pdb_file
        structure = pdb_parser.get_structure(pdb_file, os.path.join(folder,pdb_file))
        model = structure[0]
        for chain in model:
            for residue in chain:
                try:
                    ca_atoms[struct_id].append(residue['CA'])
                except:
                    ca_atoms[struct_id] = [residue['CA']]
    return ca_atoms


def find_clusters (sorted_ids,folder):
    N = len(sorted_ids)
    super_imposer = Bio.PDB.Superimposer()

    clusters_found = 0
    clusters = {clusters_found : [sorted_ids[0]]}
    
    # Read all structures CA's
    ca_atoms = get_ca_atoms(sorted_ids,folder)

    for j in sorted_ids[1:]:
        print "Peptide %d with %d.pdb" % (j, j)
        in_cluster = False
        for cluster_id in clusters.keys():
            # For each cluster representative
            representative_id = clusters[cluster_id][0]
            super_imposer.set_atoms(ca_atoms[representative_id], ca_atoms[j])
            rmsd = super_imposer.rms
            print 'RMSD between %d and %d is %5.3f' % (representative_id, j, rmsd)
            if rmsd <= 2.0:
                clusters[cluster_id].append(j)
                print "Peptide %d goes into cluster %d" % (j, cluster_id)
                in_cluster = True
                break
        
        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [j]
            print "New cluster %d" % clusters_found
    print clusters
    return clusters