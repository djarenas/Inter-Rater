#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: April 5th
"""
from inter_rater_library import *
import sys

#-----------------------------------------------------------------------------
#Input parameters
from input_parameters import *
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Read categories file; data_file and examine its properties
categories = readfile_categories(categories_filename)
array = readfile_array(data_filename)
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
rater_report(array,categories, pk_tensor) #m is number of users


#------------------------------------------------------------------------------
#Fleiss-kappa analysis. 
#------------------------------------------------------------------------------
fleiss_kappa_report(n, nij_matrix, categories)


#------------------------------------------------------------------------------
#Plot Group and Permutated kappas
#------------------------------------------------------------------------------
plot_PK_tensor(pk_tensor, nij_matrix, ymin, ymax, graph_filename, highlight, indbars)
