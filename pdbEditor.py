##Data Generating functions
def topOccupancy(PDB):
    import os, sys
    from iotbx import pdb
    from iotbx.pdb import hierarchy
    import itertools
    pdb_in = hierarchy.input(PDB)
    symm=pdb_in.crystal_symmetry()
    obj_pdb=pdb_in.construct_hierarchy()
    selected_atoms=obj_pdb.atom_selection_cache().iselection("occupancy>0.5")
    if (len(selected_atoms)>0):
        # for each in range(1,6):
        #     print (each)
        newHi=obj_pdb.select(selected_atoms)
        newHi.write_pdb_file(
            file_name=os.path.join(str("02")+"_.pdb"),
            # crystal_symmetry=symm,
            append_end=True)
    else:
        print("occupancy is lower than 0.5")



    # xray_structure=pdb_in.input.xray_structure_simple()
    # sel_cache=pdb_in.hierarchy.atom_selection_cache()
    # #selected_atoms=sel_cache.selection("occupancy<0.5")
    # print selected_atoms

        #   if (atom_group.atoms_size() != 0) :
        #     residue_group.remove_atom_group(atom_group)
        # if (residue_group.atom_groups_size() == 0) :
        #   chain.remove_residue_group(residue_group)
    # pdb_in=hierarchy.input(file_name=PDB)
    # for chains in pdb_in.hierarchy.only_model().chains():
    #     for residue_group in chains.residue_groups():
    #         for atom_group in residue_group.atom_groups():
    #             for atom in atom_group.atoms():
    #                 print atom


#highOcc=pdb_in.occupancy(>0.5)
