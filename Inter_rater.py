#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: April 5th
"""
import inter_rater_library as irl

#-----------------------------------------------------------------------------
#Input parameters
from input_file import * 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Read categories file; data_file and examine its properties
global categories
categories = irl.readfile_categories(categories_filename)
array = irl.readfile_array(data_filename)
n = len(array[0])   #number of raters
N = len(array)      #Number of subjects
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Calculate nij matrix
nij_matrix = irl.calc_nijmatrix(array,categories) 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Permutated-kappa analysis. 
#------------------------------------------------------------------------------
pk_tensor = irl.calc_permutated_kappa_tensor(array,categories)


#------------------------------------------------------------------------------
#User Statistics 
#------------------------------------------------------------------------------
irl.rater_report(array,categories, pk_tensor) #m is number of users


#------------------------------------------------------------------------------
#Fleiss-kappa analysis. 
#------------------------------------------------------------------------------
irl.fleiss_kappa_report(n, nij_matrix, categories)


#------------------------------------------------------------------------------
#Plot Group and Permutated kappas
#------------------------------------------------------------------------------
irl.plot_PK_tensor(pk_tensor, nij_matrix, ymin, ymax, graph_filename, highlight, indbars)
