#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 17:45:53 2018

@author: daniel
"""

#Data file (Rows are subjects, columns are users)
data_filename = "data_SR1.txt"

#Data file defines the ratings a rater could pick from
categories_filename = "categories.txt"

#Y-axis minimum and maximum for plot
ymin = 0
ymax = 1

#Data file name for the output graph
graph_filename = "output_figure.jpg"

#Choose which rater pairs you would like to be highlighted
highlight = ["none"]

#Option to plot error bars for each rater-pair kappas
indbars = "no"