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
    sortOccupancy("subStruc.pdb")
