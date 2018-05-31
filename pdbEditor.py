##Data Generating functions
##the following atoms2pdb function was written by Rob
def atoms2pdb(atoms):
    '''
    This function was written by Robert Offner.  It constructs a PDB hierarchy
    from a given list of atoms
    '''
    import os, sys
    from iotbx import pdb
    from iotbx.pdb import hierarchy
    r = pdb.hierarchy.root()
    r.append_model(pdb.hierarchy.model())
    r.models()[0].append_chain(pdb.hierarchy.chain(id="A"))
    #r.models()[0].chains()[0].append_residue_group(pdb.hierarchy.residue_group())
    #r.models()[0].chains()[0].residue_groups()[0].append_atom_group(pdb.hierarchy.atom_group())
    for (i,atm) in enumerate(atoms, start = 1):
        rg = pdb.hierarchy.residue_group()
        rg.resseq = i
        r.models()[0].chains()[0].append_residue_group( rg )
        ag = pdb.hierarchy.atom_group()
        rg.append_atom_group( ag )
        ag.append_atom( atm.detached_copy() )
    r.atoms_reset_serial() # make consecutive serial numbers
    return r

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

def sortOccupancy(PDB):
    '''
    This function generates PDB file of atoms having better than 0.5 occupancy.
    It sorts the atoms according to the decending order of their occupancy and
    writes a pair atoms from combination of top 5 atoms into separate PDB files.
    '''
    import os, sys
    from iotbx import pdb
    from iotbx.pdb import hierarchy
    import itertools
    mylist=[]
    pdb_in = hierarchy.input(PDB)
    symm=pdb_in.crystal_symmetry()
    obj_pdb=pdb_in.construct_hierarchy()
    selected_atoms=obj_pdb.atom_selection_cache().iselection("occupancy>0.5")
    if (len(selected_atoms)>1):
        for e in selected_atoms:
            mylist.append(obj_pdb.atoms()[e])
        sorted_atoms=sorted(mylist, key=lambda thisatom:thisatom.occ, reverse=True)
        atoms2pdb(sorted_atoms).write_pdb_file(file_name="topOcc_.pdb")
        ##the following list will generate the combination of the top 5 atoms
        iterableList=itertools.combinations(sorted_atoms[0:5],2)
        counter=int(1)
        for e in list(iterableList):
            atoms2pdb(e).write_pdb_file(file_name="combination"+str(counter)+"_.pdb",
                crystal_symmetry=pdb_in.input.crystal_symmetry(),
                append_end=True)
            counter=counter+1
    else:
        print("occupancy is lower than 0.5")
