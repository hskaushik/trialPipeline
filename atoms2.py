##o = sorted(j_alpha_hierarchy.atoms(), key= lambda thisatom:thisatom.occ, reverse=True)

def atoms2pdb(atoms):
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

