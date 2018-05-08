#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      oeffner
#
# Created:     09/01/2018
# Copyright:   (c) oeffner 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys, os
from iotbx import pdb
from cctbx import uctbx
from cctbx import sgtbx
from phaser.utils2 import *

chainid = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"



def DiffVec(vec1,vec2):
  d=[0,0,0]
  d[0] = vec1[0] - vec2[0]
  d[1] = vec1[1] - vec2[1]
  d[2] = vec1[2] - vec2[2]
  return d


def ExpandASUToP1(pdbfname, na, nb, nc, incr):
  """
  Expand a PDB crystal from the ASU to a P1 supercell with a,b and c axis
  repeating themselves ua, ub and uc times and scaled isotropically by
  1+incr along each axis.
  """
  pdbxtal = pdb.input(file_name= pdbfname)
  xtalsym = pdbxtal.crystal_symmetry()
  sg = xtalsym.space_group()
  ucell = xtalsym.unit_cell()
  xtalucell = xtalsym.unit_cell().parameters()

  bigucell = uctbx.unit_cell(( xtalucell[0]*(1.0+incr), xtalucell[1]*(1.0+incr),
   xtalucell[2]*(1.0+incr), xtalucell[3], xtalucell[4], xtalucell[5] ))
  bigxtalsym = xtalsym
  bigxtalsym._unit_cell = bigucell

  symops = sg.all_ops()
  r = pdb.hierarchy.root()
  pdb_hierarchy = pdbxtal.construct_hierarchy()
  original_sites = pdb_hierarchy.deep_copy().atoms().extract_xyz()

  P1hierar = []
  n = 0
  atoms = pdb_hierarchy.atoms()
  origsites_frac = ucell.fractionalize(sites_cart=atoms.extract_xyz())

  for a in range(na):
    for b in range(nb):
      for c in range(nc):
        newcell_sites = origsites_frac + (a, b, c)

        for (symcount, symop) in enumerate(symops):
          rot = symop.r()
          trans = symop.t()
          print "symop:(%s), ucell:[%d,%d,%d]" %(symop.as_xyz(), a, b, c )

          # put all xtal atoms back to original position
          atoms.set_xyz(newcell_sites)
          # cast back to cartesian coordinates
          atoms.set_xyz(ucell.orthogonalize(sites_frac= newcell_sites))
          # get us into fractional coordinates
          sites_frac = ucell.fractionalize(sites_cart=atoms.extract_xyz())
          # get the symmetry mate w.r.t. this space group
          new_sites = rot.as_double() * sites_frac + trans.as_double()
          # cast new origin back to unit cell around the test atoms
          origo = (0.5 + a, 0.5 + b, 0.5 + c)
          while DiffVec(new_sites.mean(), origo)[0] > 0.5:
            new_sites = new_sites - (1, 0, 0)
          while DiffVec(new_sites.mean(), origo)[0] < -0.5:
            new_sites = new_sites + (1, 0, 0)
          while DiffVec(new_sites.mean(), origo)[1] > 0.5:
            new_sites = new_sites - (0, 1, 0)
          while DiffVec(new_sites.mean(), origo)[1] < -0.5:
            new_sites = new_sites + (0, 1, 0)
          while DiffVec(new_sites.mean(), origo)[2] > 0.5:
            new_sites = new_sites - (0, 0, 1)
          while DiffVec(new_sites.mean(), origo)[2] < -0.5:
            new_sites = new_sites + (0, 0, 1)

          atoms.set_xyz(new_sites)
          # cast back to cartesian coordinates
          atoms.set_xyz(bigucell.orthogonalize(sites_frac= new_sites))
          m = pdb_hierarchy.deep_copy()
          for ch in m.chains():
            ndiv1 = n//len(chainid)
            ndiv2 = ndiv1//len(chainid)
            ch.id = chainid[n % len(chainid)]
            n += 1
          P1hierar.append(m)

  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  bigxtalucell = bigxtalsym.unit_cell().parameters()
  superucell = uctbx.unit_cell(( bigxtalucell[0]*na, bigxtalucell[1]*nb, bigxtalucell[2]*nc,
               bigxtalucell[3], bigxtalucell[4], bigxtalucell[5] ))
  superxtalsym = bigxtalsym
  superxtalsym._unit_cell = superucell
  superxtalsym._space_group_info= sgtbx.space_group_info(group=sgtbx.space_group("P 1"))

  jointpdbs = pdb.hierarchy.join_roots( P1hierar, None)
  outfname = os.path.split(pdbfname)[1].split(".pdb")[0] 
  outfname = outfname.split(".cif")[0] \
   + "[" + str(na) + "," + str(nb) + "," + str(nc) + "]_infl" + str(Roundoff(incr,6)) +".pdb"
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  #jointpdbs.write_mmcif_file(file_name=outfname, crystal_symmetry=bigxtalsym)

  refpdbout = open(outfname,"w")
  refpdbout.write(jointpdbs.as_pdb_string(bigxtalsym))
  refpdbout.close()

  return outfname



if __name__ == '__main__':
  dummy, pdbfname, na, nb, nc, incr = tuple(sys.argv)
  ExpandASUToP1(pdbfname, int(na), int(nb), int(nc), float(incr))

