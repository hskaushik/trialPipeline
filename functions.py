##Data Generating functions
def runAutosol(mtz, seq, subStruc, atom_type, wave_length,path):
    from phenix.command_line import autosol
    print("the current path is: "+path)
    if subStruc is None:
        autosol.run_autosol ([mtz, seq, atom_type, "lambda=%f" %wave_length])
    else:
        obj=autosol.run_autosol ([mtz, seq, subStruc, atom_type, "lambda=%f" %wave_length])

##SAD phaser
def runSADphaser(mtz,pdbid,PDB,seq, atom_type, wave_length,path) :
    from phaser import *
    from cctbx import xray
    print("the current path is: "+path)
    i = InputEP_DAT()
    HKLIN = mtz
    xtalid = pdbid
    waveid = "cuka"
    lamda=wave_length
    i.setHKLI(HKLIN)
    i.addCRYS_ANOM_LABI(xtalid,waveid,"F(+)","SIGF(+)","F(-)","SIGF(-)")
    i.setMUTE(False)
    r = runEP_DAT(i)
    if r.Success():
        hkl = r.getMiller()
        Fpos = r.getFpos(xtalid,waveid)
        Spos = r.getSIGFpos(xtalid,waveid)
        Ppos = r.getPpos(xtalid,waveid)
        Fneg = r.getFneg(xtalid,waveid)
        Sneg = r.getSIGFneg(xtalid,waveid)
        Pneg = r.getPneg(xtalid,waveid)
        i = InputEP_AUTO()
        i.setSPAC_HALL(r.getSpaceGroupHall())
        i.setWAVE(lamda)
        i.setCELL6(r.getUnitCell())
        i.setCRYS_MILLER(hkl)
        i.addCRYS_ANOM_DATA(xtalid,waveid,Fpos,Spos,Ppos,Fneg,Sneg,Pneg)
        i.setATOM_PDB(xtalid,PDB)
        i.setLLGC_COMP(True)
        i.addLLGC_SCAT(atom_type)
        i.addCOMP_PROT_SEQ_NUM(seq,1.)
        i.setTITL(xtalid+path)
        i.setROOT(xtalid+path)
        r = runEP_AUTO(i)
        f=open(path+'.log','w')
        f.write(r.logfile())
        f.close()
        print "LogLikelihood = " , r.getLogLikelihood()
    else:
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
    return;

##this runs the MR_SAD
def runMRsad(mtz,pdbid,PDB,seq,atom_type,wave_length,path) :
    from cStringIO import StringIO
    import pickle
    from phaser import *
    from cctbx import xray
    import os, sys
    print("the current path is: "+path)
    # original=sys.stdout
    # sys.stdout=open('log.txt','w')
    i = InputEP_DAT()
    HKLIN = mtz
    xtalid = pdbid
    lamda=wave_length
    waveid = "cuka"
    i.setHKLI(HKLIN)
    i.addCRYS_ANOM_LABI(xtalid,waveid,"F(+)","SIGF(+)","F(-)","SIGF(-)")
    i.setMUTE(True)
    # o.setPackagePhenix(file_object=redirect_str)
    r = runEP_DAT(i)
    # r = runEP_DAT(i)
    # print (r.logfile())
    if r.Success():
        hkl = r.getMiller()
        Fpos = r.getFpos(xtalid,waveid)
        Spos = r.getSIGFpos(xtalid,waveid)
        Ppos = r.getPpos(xtalid,waveid)
        Fneg = r.getFneg(xtalid,waveid)
        Sneg = r.getSIGFneg(xtalid,waveid)
        Pneg = r.getPneg(xtalid,waveid)
        i = InputEP_AUTO()
        i.setMUTE(True)
        i.setSPAC_HALL(r.getSpaceGroupHall())
        i.setWAVE(lamda)
        i.setCELL6(r.getUnitCell())
        i.setCRYS_MILLER(hkl)
        i.addCRYS_ANOM_DATA(xtalid,waveid,Fpos,Spos,Ppos,Fneg,Sneg,Pneg)
        i.setPART_PDB(PDB)
        i.setPART_DEVI(0.8)
        i.setPART_VARI("RMS")
        i.setLLGC_COMP(True)
        i.addLLGC_SCAT(atom_type)
        i.addCOMP_PROT_SEQ_NUM(seq,1.)
        i.setTITL(xtalid)
        i.setROOT(xtalid+path)
        r = runEP_SAD(i)
        print "LogLikelihood = " , r.getLogLikelihood()
        f=open(path+'.log','w')
        f.write(r.logfile())
        f.close()


        ##load r object into pickle file
        # with open("result.pkl", "w") as result:
        #     pickle.dump(r,result)
        #use the following to retrieve the pickle
        #result=pickle.load(open("result.pkl", "rb"))
    else:
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
    # sys.stdout=original
    return;
    #this can be run by following command
    #runMRsad("/home/ksh40/work/lysozyme/3abc/anomalous/french_wilson.mtz","3abc","/home/ksh40/work/lysozyme/lyso_SSAD/refined18_lysozyme.pdb","/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa",1.54)

##this code will run Phassade
def runPhassade(mtz,pdbid,seq, atom_type, wave_length,path) :
    from phaser import *
    from cctbx import xray
    path=path
    print("the current path is: "+path)
    i = InputEP_DAT()
    HKLIN = mtz
    xtalid = pdbid
    waveid = "cuka"
    lamda=wave_length
    i.setHKLI(HKLIN)
    i.addCRYS_ANOM_LABI(xtalid,waveid,"F(+)","SIGF(+)","F(-)","SIGF(-)")
    i.setMUTE(False)
    r = runEP_DAT(i)
    if r.Success():
        hkl = r.getMiller()
        Fpos = r.getFpos(xtalid,waveid)
        Spos = r.getSIGFpos(xtalid,waveid)
        Ppos = r.getPpos(xtalid,waveid)
        Fneg = r.getFneg(xtalid,waveid)
        Sneg = r.getSIGFneg(xtalid,waveid)
        Pneg = r.getPneg(xtalid,waveid)
        i = InputEP_SSD()
        i.setMUTE(True)
        i.setSPAC_HALL(r.getSpaceGroupHall())
        i.setWAVE(lamda)
        i.setCELL6(r.getUnitCell())
        i.setCRYS_MILLER(hkl)
        i.addCRYS_ANOM_DATA(xtalid,waveid,Fpos,Spos,Ppos,Fneg,Sneg,Pneg)
        i.setLLGC_COMP(True)
        i.addLLGC_SCAT(atom_type)
        i.addCOMP_PROT_SEQ_NUM(seq,1.)
        i.setTITL(xtalid)
        i.setROOT(xtalid+path)
        r = runEP_SSD(i)
        f=open(path+'.log','w')
        f.write(r.logfile())
        f.close()
    else:
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
    return;

##this code will run Shelx pipeline
def runShelx(mtz,pdbid,seq,atom_type,wave_length,solvent_content,resolution,path):
    import os
    assert os.path.isfile(mtz)
    from iotbx import crystal_symmetry_from_any
    print("the current path is: "+path)
    crystal_data=crystal_symmetry_from_any.extract_from(mtz)
    solvent_content=str(solvent_content)
    resolution=str(resolution)
    unit_cell=str(crystal_data.unit_cell()).replace(",","")
    unit_cell=unit_cell.replace("(","")
    unit_cell=unit_cell.replace(")","")
    spacegroup=str(crystal_data.space_group_info()).replace(" ","")
    os.system('mtz2sca ' +mtz)
    sca= pdbid+'.sca'
    pdbid=pdbid+path
    copyFunction= 'cp /home/ksh40/work/trial_runs/sampleScripts/shelx_sad.sh '+pdbid+ '.sh'
    os.system(copyFunction)
    #replace regular expressions
    sedReplace='sed -i \'s/replace_pdbid/'+pdbid+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplaceFa='sed -i \'s/replace_pdbid_fa/'+pdbid+'_fa/\' ' +pdbid+'.sh'
    os.system(sedReplaceFa)
    sedReplace='sed -i \'s/replace_solvent_content/'+solvent_content+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplace='sed -i \'s/replace_resolution/'+resolution+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplace='sed -i \'s/replace_dimension/'+unit_cell+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplace='sed -i \'s/replace_spacegroup/'+spacegroup+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplace='sed -i \'s/replace_atomtype/'+atom_type+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    sedReplace='sed -i \'s/replace_input_data/'+sca+'/\' ' +pdbid+'.sh'
    os.system(sedReplace)
    os.system('sh ' +pdbid+'.sh' + ">"+path+"shelx.log" )

##Analysing functions
def runEmma(referencePDB, currentPDB,path):
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
    from phenix.command_line import get_cc_mtz_pdb
    print("this is currentMTZ: ", currentMTZ)
    print("this is referencePDB: ", referencePDB)
    obj=get_cc_mtz_pdb.get_cc_mtz_pdb([currentMTZ, referencePDB])

def runExpand2P1(PDB):
    import ExpandASU
    ExpandASU.ExpandASUToP1(PDB, 1, 1, 1, 0)
