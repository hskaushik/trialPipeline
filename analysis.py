##Analysing functions
def runEmma(referencePDB, currentPDB,path):
    '''
    This function compares two substructures
    '''
    from iotbx import crystal_symmetry_from_any
    import iotbx.pdb
    from iotbx.cns import sdb_reader
    from iotbx.kriber import strudat
    from iotbx.option_parser import option_parser
    from cctbx import euclidean_model_matching as emma
    import sys, os
    import cctbx.xray
    from iotbx.command_line import emma
    emma.run(referencePDB, currentPDB)

def runCC_MTZ_PDB(currentMTZ, referencePDB,path):
    '''
    This function generates correlation coefficient by comparing structure
    factors against a PDB model
    '''
    from phenix.command_line import get_cc_mtz_pdb
    print("this is currentMTZ: ", currentMTZ)
    print("this is referencePDB: ", referencePDB)
    obj=get_cc_mtz_pdb.get_cc_mtz_pdb([currentMTZ, referencePDB])

def runExpand2P1(PDB):
    '''
    This function was provided by Rober Offner.  It expands the structure into
    P1 space group.
    '''
    import ExpandASU
    ExpandASU.ExpandASUToP1(PDB, 1, 1, 1, 0)

# def atom_count(PDB):
#     '''
#     this function will count the number of atoms in a given PDB file
#     '''
#     from iotbx.pdb import hierarchy
#     pdb_in=hierarchy.input(PDB)
#     obj_pdb=pdb_in.construct_hierarchy()
#     selected_atoms=obj_pdb.atom_selection_cache().iselection("occupancy>0.1")
#     print (len(selected_atoms))
