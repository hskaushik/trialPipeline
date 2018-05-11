from __future__ import division
import sys
import os
import csv
import glob
from functions import *
from pdbEditor import sortOccupancy

# class Job_Info:
#     pdbid=None


if (__name__ == "__main__"):
    # topOccupancy("/home/ksh40/work/trial_runs/2od6/runs/pc/SADphaser/2od6_pc_mrSAD.1.pdb")
    # sortOccupancy("autosol.pdb")

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

    dirs = (
        ["runs"],
        ["runs", "pc"],
        ["runs", "pc","mrSAD"],
        ["runs", "pc","mrSAD","SADphaser"],
        ["runs", "pc","SADphaser"],
        ["runs", "pc","autosol"],
        ["runs", "default"],
        ["runs", "default","shelx"],
        ["runs", "default","shelx","SADphaser"],
        ["runs", "default","autosol"],
        ["runs", "default","phassade"],
        ["runs", "default","phassade","SADphaser"],
        )

    for dir_list in dirs:
        pathname = os.path.join(home_path, *dir_list)
        if not os.path.exists(pathname):
            os.makedirs(pathname)

    # if not os.path.exists(os.path.join(home_path,"runs")):
    #     os.makedirs("runs")

    # ##creating directory hierarchy
    # if not os.path.exists(os.path.join(home_path,"runs","pc")):
    #     os.makedirs(os.path.join(home_path,"runs","pc"))
    # if not os.path.exists(os.path.join(home_path,"runs","pc","mrSAD")):
    #     os.makedirs(os.path.join(home_path,"runs","pc","mrSAD"))
    # if not os.path.exists(os.path.join(home_path,"runs","pc","SADphaser")):
    #     os.makedirs(os.path.join(home_path,"runs","pc","SADphaser"))
    # if not os.path.exists(os.path.join(home_path,"runs","pc","autosol")):
    #     os.makedirs(os.path.join(home_path,"runs","pc","autosol"))
    # if not os.path.exists(os.path.join(home_path,"runs","default")):
    #     os.makedirs(os.path.join(home_path,"runs","default"))
    # if not os.path.exists(os.path.join(home_path,"runs","default","autosol")):
    #     os.makedirs(os.path.join(home_path,"runs","default","autosol"))
    # if not os.path.exists(os.path.join(home_path,"runs","default","phassade")):
    #     os.makedirs(os.path.join(home_path,"runs","default","phassade"))
    # if not os.path.exists(os.path.join(home_path,"runs","default","phassade","SADphaser")):
    #     os.makedirs(os.path.join(home_path,"runs","default","phassade","SADphaser"))
    # if not os.path.exists(os.path.join(home_path,"runs","default","shelx")):
    #     os.makedirs(os.path.join(home_path,"runs","default","shelx"))
    # if not os.path.exists(os.path.join(home_path,"runs","default","shelx","SADphaser")):
    #     os.makedirs(os.path.join(home_path,"runs","default","shelx","SADphaser"))


    # ##starting positive control pipeline
    os.chdir(os.path.join(home_path,"runs","pc","mrSAD"))
    path="_pc"+"_mrSAD"
    # runMRsad(mtz,pdbid,PDB,seq,"AX",resolution,path)
    mrSAD_subStruc=os.path.join(home_path,"runs","pc","mrSAD",pdbid+path+".1.pdb")

    ### spliting the substructure into combination of multiple PDB files
    sortOccupancy(mrSAD_subStruc)

    ##starting SADphaser inside mrSAD
    os.chdir(os.path.join(home_path,"runs","pc","mrSAD","SADphaser"))
    path="_pc"+"_mrSAD"+"_SADphaser"
    # runSADphaser(mtz,pdbid,mrSAD_subStruc,seq, atom_type, wave_length, path)
    ##starting multiple SADphaser inside mrSAD for subset of substructure
    os.chdir(os.path.join(home_path,"runs","pc","mrSAD"))
    current_dir=os.getcwd()
    for file_path in glob.glob(os.path.join(current_dir,'*_.pdb')):
        new_dir=file_path.rsplit('.',1)[0]
        os.mkdir(os.path.join(new_dir))
        os.chdir(os.path.join(new_dir))
        print ("this is new_dir %s", new_dir)
        path="_pc"+"_mrSAD"+"_SADphaser"+"_"+str(new_dir.split("/")[-1])
        print ("this is path %s", path)
        mrSAD_subStruc=file_path
        print("this is file_path", file_path)
        runSADphaser(mtz,pdbid,mrSAD_subStruc,seq, atom_type, wave_length, path)
        os.chdir(os.path.join(current_dir))



    # ##starting positive control phaser
    os.chdir(os.path.join(home_path,"runs","pc","SADphaser"))
    original_subStruc=os.path.join(home_path,"subStruc.pdb")
    path="_pc"+"_SADphaser"
    # runSADphaser(mtz,pdbid,original_subStruc,seq, atom_type, wave_length, path)
    ### spliting the substructure into combination of multiple PDB files
    sortOccupancy(original_subStruc)

    ##starting multiple SADphaser inside using original substructure
    # os.chdir(os.path.join(home_path,"runs","pc","SADphaser"))
    current_dir=os.getcwd()
    for file_path in glob.glob(os.path.join(current_dir,'*_.pdb')):
        new_dir=file_path.rsplit('.',1)[0]
        os.mkdir(os.path.join(new_dir))
        os.chdir(os.path.join(new_dir))
        print ("this is new_dir %s", new_dir)
        path="_pc"+"_SADphaser"+"_"+str(new_dir.split("/")[-1])
        print ("this is path %s", path)
        original_subStruc=file_path
        print("this is file_path", file_path)
        runSADphaser(mtz,pdbid,original_subStruc,seq, atom_type, wave_length, path)
        os.chdir(os.path.join(current_dir))

    #
    # ##starting Autosol for positive control
    # os.chdir(os.path.join(home_path,"runs","pc","autosol"))
    # path="_pc"+"_Autosol"
    # runAutosol(mtz, seq, original_subStruc, atom_type, wave_length,path)
    #
    #
    # ##starting default pipeline
    # ##starting autosol with no substructure information
    # os.chdir(os.path.join(home_path,"runs","default","autosol"))
    # path="_default"+"_Autosol"
    # runAutosol(mtz, seq, None, atom_type, wave_length,path)
    #
    # ##starting shelx under default pipeline
    os.chdir(os.path.join(home_path,"runs","default","shelx"))
    path="_default"+"_shelx"
    runShelx(mtz,pdbid,seq,atom_type,wave_length,solvent_content,resolution,path)
    shelx_subStruc=os.path.join(home_path,"runs","default","shelx",pdbid+path+"_fa.pdb")

    ### spliting the substructure into combination of multiple PDB files
    sortOccupancy(shelx_subStruc)

    ##starting multiple SADphaser inside using original substructure
    current_dir=os.getcwd()
    for file_path in glob.glob(os.path.join(current_dir,'*_.pdb')):
        new_dir=file_path.rsplit('.',1)[0]
        os.mkdir(os.path.join(new_dir))
        os.chdir(os.path.join(new_dir))
        print ("this is new_dir %s", new_dir)
        path="_default"+"_shelx"+"_"+str(new_dir.split("/")[-1])
        print ("this is path %s", path)
        shelx_subStruc=file_path
        print("this is file_path", file_path)
        runSADphaser(mtz,pdbid,shelx_subStruc,seq, atom_type, wave_length, path)
        os.chdir(os.path.join(current_dir))

    #
    # ##starting SADphaser inside shelx
    # os.chdir(os.path.join(home_path,"runs","default","shelx","SADphaser"))
    # shelx_subStruc=os.path.join("../",pdbid+path+"_fa.pdb")
    # path="_default"+"_shelx"+"_SADphaser"
    # runSADphaser(mtz,pdbid,shelx_subStruc,seq, atom_type, wave_length, path)

# shelxdatas = []
# for f in filenames:
#     shelxdatas.append(ShelxData(filename))




# class ShelxData:
#     '''
#     This class does...
#     '''
#     def __init__(self, filename):
#         '''
#         This function does...
#         '''
#         self._parse_file()
#
#
#     def _parse_file(self, filename):
#         pass
#
