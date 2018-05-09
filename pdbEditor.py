##Data Generating functions
def topOccupancy(PDB):
    import os, sys
    from iotbx import pdb
    from iotbx.pdb import hierarchy
    import itertools
    occ=float(0.1)
    pdb_in = hierarchy.input(PDB)
    symm=pdb_in.crystal_symmetry()
    obj_pdb=pdb_in.construct_hierarchy()
    selected_atoms=obj_pdb.atom_selection_cache().iselection("occupancy>"+str(occ)+" ")
    counter=int(len(selected_atoms))
    # counter=17
    if (counter>0):
        while (counter>2):
            occ=occ+float(0.02)
            # print ("value of counter is %d and occ is %f",counter, occ)
            selected_atoms=obj_pdb.atom_selection_cache().iselection("occupancy>"+str(occ)+" ")
            counter=int(len(selected_atoms))
            if (counter<6):
                # print ("value of counter inside if is %d ",counter)
                newHi=obj_pdb.select(selected_atoms)
                newHi.write_pdb_file(file_name=os.path.join(str(counter)+"_.pdb"))

    else:
        print("occupancy is lower than 0.1")
