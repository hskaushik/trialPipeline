from __future__ import division
import sys
import os
import csv
from functions import *
# from pdbEditor import topOccupancy

if (__name__ == "__main__"):
    # topOccupancy("/home/ksh40/work/trial_runs/2od6/runs/pc/SADphaser/2od6_pc_mrSAD.1.pdb")

    # runMRsad("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","2-7","/home/ksh40/work/dorothee/pdbDump/5wzq.pdb","/home/ksh40/work/dorothee/pdbDump/5wzq.fa","AX",2.7, "trial")
    # runSADphaser("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","2-7","2-7trial.1.pdb","/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","S", 2.7, "trialSAD")
    # runShelx("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","NagBb_2-7", "/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","s",2.7,"0.5",2.4,"trialShelx")
    # runShelx("/home/ksh40/work/dorothee/NagBb_2-7.mtz","NagBb_2-7","/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","s",2.7,"0.5",2.4,"")
    # runPhassade("/home/ksh40/work/trial_runs/2od6/runs/trial/NagBb_2-7.mtz","NagBb_2-7","/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","s", 2.7,"trialPhassade")

    home_path=os.getcwd()

    ##read from the csv file
    with open('input.data', 'r') as csvfile:
        read_list=csv.DictReader(csvfile)
        for row in read_list:

            ##define standard variables
            pdbid=row['pdbid']
            atom_type=row['atom_type']
            wave_length=float(row['wave_length'])
            solvent_content=float(row['solvent_content'])
            resolution=float(row['resolution'])

            ##define variables with paths
            PDB=os.path.join(home_path,pdbid+".pdb")
            seq=os.path.join(home_path,pdbid+".fa")
            mtz=os.path.join(home_path,pdbid+".mtz")
            print("path to PDB model file: "+PDB)
            print("path to sequence fasta file : "+seq)
            print("path to FW pre-converted MTZ file : "+mtz)

    if not os.path.exists(os.path.join(home_path,"runs")):
        os.makedirs("runs")

    ##creating directory hierarchy
    if not os.path.exists(os.path.join(home_path,"runs","pc")):
        os.makedirs(os.path.join(home_path,"runs","pc"))
    if not os.path.exists(os.path.join(home_path,"runs","pc","mrSAD")):
        os.makedirs(os.path.join(home_path,"runs","pc","mrSAD"))
    if not os.path.exists(os.path.join(home_path,"runs","pc","SADphaser")):
        os.makedirs(os.path.join(home_path,"runs","pc","SADphaser"))
    if not os.path.exists(os.path.join(home_path,"runs","pc","autosol")):
        os.makedirs(os.path.join(home_path,"runs","pc","autosol"))
    if not os.path.exists(os.path.join(home_path,"runs","default")):
        os.makedirs(os.path.join(home_path,"runs","default"))
    if not os.path.exists(os.path.join(home_path,"runs","default","autosol")):
        os.makedirs(os.path.join(home_path,"runs","default","autosol"))
    if not os.path.exists(os.path.join(home_path,"runs","default","phassade")):
        os.makedirs(os.path.join(home_path,"runs","default","phassade"))
    if not os.path.exists(os.path.join(home_path,"runs","default","phassade","SADphaser")):
        os.makedirs(os.path.join(home_path,"runs","default","phassade","SADphaser"))
    if not os.path.exists(os.path.join(home_path,"runs","default","shelx")):
        os.makedirs(os.path.join(home_path,"runs","default","shelx"))
    if not os.path.exists(os.path.join(home_path,"runs","default","shelx","SADphaser")):
        os.makedirs(os.path.join(home_path,"runs","default","shelx","SADphaser"))


    ##starting positive control pipeline
    os.chdir(os.path.join(home_path,"runs","pc","mrSAD"))
    path="_pc"+"_mrSAD"
    runMRsad(mtz,pdbid,PDB,seq,"AX",resolution,path)

    ##starting SADphaser inside mrSAD
    if not os.path.exists(os.path.join(home_path,"runs","pc","mrSAD","SADphaser")):
        os.makedirs(os.path.join(home_path,"runs","pc","mrSAD","SADphaser"))
    os.chdir(os.path.join(home_path,"runs","pc","mrSAD","SADphaser"))
    mrSAD_subStruc=os.path.join("../",pdbid+path+".1.pdb")
    path="_pc"+"_mrSAD"+"_SADphaser"
    runSADphaser(mtz,pdbid,mrSAD_subStruc,seq, atom_type, wave_length, path)

    ##starting positive control phaser
    os.chdir(os.path.join(home_path,"runs","pc","SADphaser"))
    original_subStruc=os.path.join(home_path,"subStruc.pdb")
    path="_pc"+"_SADphaser"
    runSADphaser(mtz,pdbid,original_subStruc,seq, atom_type, wave_length, path)

    ##starting Autosol for positive control
    os.chdir(os.path.join(home_path,"runs","pc","autosol"))
    path="_pc"+"_Autosol"
    runAutosol(mtz, seq, original_subStruc, atom_type, wave_length,path)


    ##starting default pipeline
    ##starting autosol with no substructure information
    os.chdir(os.path.join(home_path,"runs","default","autosol"))
    path="_default"+"_Autosol"
    runAutosol(mtz, seq, None, atom_type, wave_length,path)

    ##starting shelx under default pipeline
    os.chdir(os.path.join(home_path,"runs","default","shelx"))
    path="_default"+"_shelx"
    runShelx(mtz,pdbid,seq,atom_type,wave_length,solvent_content,resolution,path)

    ##starting SADphaser inside shelx
    if not os.path.exists(os.path.join(home_path,"runs","default","shelx","SADphaser")):
        os.makedirs(os.path.join(home_path,"runs","default","shelx","SADphaser"))
    os.chdir(os.path.join(home_path,"runs","default","shelx","SADphaser"))
    shelx_subStruc=os.path.join("../",pdbid+path+"_fa.pdb")
    path="_default"+"_shelx"+"_SADphaser"
    runSADphaser(mtz,pdbid,shelx_subStruc,seq, atom_type, wave_length, path)



    # runMRsad("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","2-7","/home/ksh40/work/dorothee/pdbDump/5wzq.pdb","/home/ksh40/work/dorothee/pdbDump/5wzq.fa","AX",2.7, "trial")
    # #runSADphaser("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","2-7","../2-7.1.pdb","/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","S", 2.7)
    # #runAutosol("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz", "/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa", "/home/ksh40/work/dorothee/pdbDump/subStruc.pdb", "S", 2.7)
    # #runShelx("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz","NagBb_2-7", "/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa","s",2.7,"0.5",2.4)
    # #runnning autosol as default without substructure information
    # #runAutosol("/home/ksh40/work/dorothee/NagBb_2-7_fw.mtz", "/home/ksh40/work/lysozyme/lyso_SSAD/seq.fa", None, "S", 2.7)
    # #runCC_MTZ_PDB("2-7.1.mtz", "final.pdb")
