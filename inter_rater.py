#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: October2022
"""
import inter_rater_library as irl

#-----------------------------------------------------------------------------
#Input parameters
from input_file import * 
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Read categories file and data_file; check their properties
#------------------------------------------------------------------------------
print("Reading and checking data from files...")
categories = irl.read_file(categories_filename)
nparray = irl.read_file(data_filename)

#Terminate if the data is not usable for inter-rater analysis.
if (not irl.check_input_data(categories, nparray)):
    quit()
else:
    print("No errors found in the input data files.")


#------------------------------------------------------------------------------
#Print Statistics For Each Rater 
#------------------------------------------------------------------------------
permutated_kappa = irl.calc_PermutatedKappa(nparray, categories)
irl.print_RaterReports(nparray, categories, permutated_kappa)


#------------------------------------------------------------------------------
#Print Fleiss-kappa analysis. 
#------------------------------------------------------------------------------
n = len(nparray[0]) #Number of users
nij = irl.calc_nijmatrix(nparray, categories)
irl.print_GroupReport(n, nij, categories)


#------------------------------------------------------------------------------
#Plot Group and Permutated kappas
#------------------------------------------------------------------------------
kappa = irl.calc_Kappa(nij)
kvariance = irl.estimate_KappaVariance(nij)
irl.plot_Kappas((kappa, kvariance**0.5), permutated_kappa, \
                graph_filename, (ymin,ymax))