#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Inter-Rater
Version: 1.0
@author: Daniel J. Arenas
Date: 04/24/2018
Submitted to Journal of Open Research Software (JORS) as:
'Inter-rater: Software for analysis of inter-rater reliability by permutating
pairs of multiple users.' 
"""

from Inter_rater_library import *
import argparse
import sys

#------------------------------------------------------------------------------
#Arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument("-dfile", dest="data_filename", default = "missing")
parser.add_argument("-cfile", dest="categories_filename", default = "categories.txt")
parser.add_argument("-ofile", dest="graph_filename", default = "output_figure.jpg")
parser.add_argument("-ymin", type = float, dest ="ymin", default = 0)
parser.add_argument("-ymax", type = float, dest ="ymax", default = 1)
parser.add_argument("-highlight", dest = "highlight", default = ["none"])
parser.add_argument("-indbars", dest = "indbars", default = "no")

args = parser.parse_args()
if args.data_filename == "missing":
    print("You must enter a data file to be analyzed")
    sys.exit()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Read categories file; data_file and examine its properties
categories = ReadFile_Categories(args.categories_filename)
array = ReadFile_Array(args.data_filename)
n = len(array[0])   #number of raters
N = len(array)      #Number of subjects
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Calculate nij matrix
nij_matrix = Calculate_nijmatrix(array,categories) 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Permutated-kappa analysis. 
#------------------------------------------------------------------------------
pk_tensor = Calculate_permutated_kappa_tensor(array,categories)


#------------------------------------------------------------------------------
#Print out user statistics 
#------------------------------------------------------------------------------
Stats(array,categories, pk_tensor) #m is number of users


#------------------------------------------------------------------------------
#Print out group statistics. (Fleiss-kappa analysis.) 
#------------------------------------------------------------------------------
Fleiss_kappa_analysis(n, nij_matrix, categories)


#------------------------------------------------------------------------------
#Plot Group and Permutated kappas
#------------------------------------------------------------------------------
Plot_PK_tensor(pk_tensor, nij_matrix, args.ymin, args.ymax, \
               args.graph_filename, args.highlight, args.indbars)