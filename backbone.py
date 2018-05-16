from __future__ import division
import sys, os, csv
import glob
import multiprocessing
from multiprocessing import Pool
from functions import *
from pdbEditor import sortOccupancy


# class Job_Info:
#     pdbid=None

if (__name__ == "__main__"):
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
            nproc=int(10) # number of CPUs to be used for parallel processing

            ##define variables with paths
            PDB=os.path.join(home_path,pdbid+".pdb")
            seq=os.path.join(home_path,pdbid+".fa")
            mtz=os.path.join(home_path,pdbid+".mtz")
            subStruc=os.path.join(home_path,"subStruc.pdb")
            print("path to PDB model file: "+PDB)
            print("path to sequence fasta file : "+seq)
            print("path to FW pre-converted MTZ file : "+mtz)
            print("path to original substructure is : "+subStruc)

    dirs = (
        ["runs"],
        ["runs", "pc"],
        ["runs", "pc","mrSAD"],
        ["runs", "pc","SADphaser"],
        ["runs", "pc","autosol"],
        ["runs", "default"],
        ["runs", "default","shelx"],
        ["runs", "default","autosol"],
        ["runs", "default","phassade"],
        )

    for dir_list in dirs:
        full_path = os.path.join(home_path, *dir_list)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        os.chdir(full_path)
        path_name=(full_path.rsplit('runs',1)[1]).replace("/","_")
        print("path_name is", path_name)
        if path_name=="_pc_mrSAD":
            ##run mrSAD
            try:
                # runMRsad(mtz,pdbid,PDB,seq,"AX",resolution,path_name)
                mrSAD_subStruc=os.path.join(full_path,pdbid+path_name+".1.pdb")

                ### spliting the substructure into combination of multiple PDB files
                sortOccupancy(mrSAD_subStruc)

                ##run SADphaser inside mrSAD
                if not os.path.exists(os.path.join(full_path,"SADphaser")):
                    os.makedirs(os.path.join(full_path,"SADphaser"))
                os.chdir(os.path.join(full_path,"SADphaser"))
                path_name=(path_name+"_SADphaser")
                # runSADphaser(mtz,pdbid,mrSAD_subStruc,seq, atom_type, wave_length, path_name,None)
                os.chdir(full_path)
        #
        #         ##running SADphaser for each of the top combinations
        #         SADphaser_args=[]
        #         for file_name in glob.glob('*_.pdb'):
        #             subStruc_combination_file_path=os.path.join(full_path,file_name)
        #             subStruc_path_name=(path_name+"_"+str(file_name.rsplit('.',1)[0]))
        #             dir_path=os.path.join(full_path,str(file_name.rsplit('.',1)[0]))
        #             SADphaser_args.append([mtz,pdbid,subStruc_combination_file_path,seq,atom_type,wave_length,subStruc_path_name,dir_path])
        #
        #         ##run SADphaser in parallel using multiprocessing
        #         # pool=Pool(nproc)
        #         # pool.map(multi_run_wrapper,SADphaser_args)
            except:
                pass
        if path_name=="_pc_SADphaser":
            ##run SADphaser using the original substructure
            try:
                # runSADphaser(mtz,pdbid,subStruc,seq, atom_type, wave_length, path_name,None)
            except:
                pass
        if path_name=="_default_shelx":
            ##run shelx
            try:
                # runShelx(mtz,pdbid,seq,atom_type,wave_length,solvent_content,resolution,path_name)
                shelx_subStruc=os.path.join(full_path,pdbid+path_name+"_fa.pdb")

                ### spliting the substructure into combination of multiple PDB files
                sortOccupancy(shelx_subStruc)

                ##run SADphaser inside mrSAD
                if not os.path.exists(os.path.join(full_path,"SADphaser")):
                    os.makedirs(os.path.join(full_path,"SADphaser"))
                os.chdir(os.path.join(full_path,"SADphaser"))
                path_name=(path_name+"_SADphaser")
                # runSADphaser(mtz,pdbid,shelx_subStruc,seq, atom_type, wave_length, path_name, None)
                os.chdir(full_path)

                ##running SADphaser for each of the top combinations
                SADphaser_args=[]
                for file_name in glob.glob('*_.pdb'):
                    subStruc_combination_file_path=os.path.join(full_path,file_name)
                    subStruc_path_name=(path_name+"_"+str(file_name.rsplit('.',1)[0]))
                    dir_path=os.path.join(full_path,str(file_name.rsplit('.',1)[0]))
                    SADphaser_args.append([mtz,pdbid,subStruc_combination_file_path,seq,atom_type,wave_length,subStruc_path_name,dir_path])

                ##run SADphaser in parallel using multiprocessing
                pool=Pool(nproc)
                pool.map(multi_run_wrapper,SADphaser_args)
            except:
                pass
        if path_name=="_default_autosol":
            ##run autosol without substructure information
            try:
                runAutosol(mtz, seq, None, atom_type, wave_length,path_name)
            except:
                pass
        if path_name=="_default_phassade":
            ##run phassade
            try:
                runPhassade(mtz,pdbid,seq,atom_type,wave_length,path_name)
                SADphaser_args=[]
                ##for each solution obtained from phassade
                for file_name in glob.glob('*_phassade.*.pdb'):
                    subStruc_combination_file_path=os.path.join(full_path,file_name)
                    subStruc_path_name=(path_name+"_"+str(file_name.rsplit('.',1)[0]))
                    dir_path=os.path.join(full_path,str(file_name.rsplit('.',1)[0]))
                    SADphaser_args.append([mtz,pdbid,subStruc_combination_file_path,seq,atom_type,wave_length,subStruc_path_name,dir_path])

                ##run SADphaser in parallel using multiprocessing
                pool=Pool(nproc)
                pool.map(multi_run_wrapper,SADphaser_args)
            except:
                pass

# # class ShelxData:
# #     '''
# #     This class does...
# #     '''
# #     def __init__(self, filename):
# #         '''
# #         This function does...
# #         '''
# #         self._parse_file()
# #
# #
# #     def _parse_file(self, filename):
# #         pass
# #
