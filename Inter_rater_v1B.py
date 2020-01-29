#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 13:20:59 2017
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: April 9th, 2018 in Osaka
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
categories = read_categories_file(args.categories_filename)
array = read_arrays_file(args.data_filename)
n = len(array[0])   #number of raters
N = len(array)      #Number of subjects
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Calculate nij matrix
nij_matrix = calculate_nijmatrix(array,categories) 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Permutated-kappa analysis. 
#------------------------------------------------------------------------------
pk_tensor = calculate_permutated_kappa_tensor(array,categories)


#------------------------------------------------------------------------------
#User Statistics 
#------------------------------------------------------------------------------
calculate_user_stats(array,categories, pk_tensor) #m is number of users


#------------------------------------------------------------------------------
#Fleiss-kappa analysis. 
#------------------------------------------------------------------------------
fleiss_kappa_analyze(n, nij_matrix, categories)


#------------------------------------------------------------------------------
#Plot Group and Permutated kappas
#------------------------------------------------------------------------------
plot_PK_tensor(pk_tensor, nij_matrix, args.ymin, args.ymax, args.graph_filename, args.highlight, args.indbars)
