#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Title: Library of Inter-Rater 
Version: 1.0
@author: Daniel J. Arenas
Date: 04/24/2018
Submitted to Journal of Open Research Software (JORS) as:
'Inter-rater: Software for analysis of inter-rater reliability by permutating
pairs of multiple users.' 
"""
#=====================Import===================================
#==============================================================
import csv
import math
import matplotlib.pyplot as plt
#==============================================================
#==============================================================


#==================Functions for ReadingFiles==========
#==============================================================
def ReadFile_Categories(filename):
    array = []        #Lists are mutable
    with open(filename, 'r') as inF: 
      for line in csv.reader(inF, dialect="excel-tab"):
            array = array + line
    return(array)

def ReadFile_Array(filename):
    array = []        #Lists are mutable
    with open(filename, 'r') as inF:
        for line in csv.reader(inF, dialect="excel-tab"):
            array = array + [line]
    return(array)
#==============================================================
#==============================================================


#====================Functions for nij matrix==================
#==============================================================
def Calculate_nij(array, i, rating):
    """
    i: study number  (1 to N)
    """
    n = len(array[0])   #Number of raters
    row = array[i]
    nij = 0
    for l in range(0,n):
        if row[l] == rating:
            nij = nij + 1
    return(nij)

def Calculate_row(array,i,categories):
    """
    i: study number (1 to N)
    categories: list of possible categories "-1,0,1"
    """
    k = len(categories)   #Number of categories
    row = []
    for m in range(0,k):
        nij = Calculate_nij(array,i,categories[m])
        row = row + [nij]
    return(row)

def Calculate_nijmatrix(array,categories):
    """
    array format:
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1
    output: nij matrix (See Appendix A of JORS publication)
    """
    matrix = []
    N = len(array)
    for m in range(0,N):
        row = Calculate_row(array,m,categories)
        matrix = matrix + [row]
    return(matrix)
#==============================================================
#==============================================================


#===============Functions for Pi column and pj row=============
#==============================================================
def Calculate_Pi(matrix, i):
    """ 
    nij matrix (See Appendix A of software's JORS publication)
    i: study number
    Calculate agreement for study i
    """
    k = len(matrix[i])
    n = 0.0
    for m in range(0,k):
        n = n + matrix[i][m]
    total = 0.0
    for m in range(0,k):
        nij = matrix[i][m]
        total = total + (1/(n**2-n))*(nij**2-nij)
    return(total)

def Calculate_Pi_column(matrix):
    N = len(matrix)    
    Pi_column = []
    for l in range(0,N):
        value = Calculate_Pi(matrix, l)
        Pi_column = Pi_column + [value]
    return(Pi_column)

def Calculate_P(matrix):
    Pi_column = Calculate_Pi_column(matrix)
    return(sum(Pi_column)/float(len(matrix)))

def Calculate_Total_row(matrix):
    Total_row = []    
    for i in range (0, len(matrix[0])):
        total = 0
        for m in range (0,len(matrix)):
            total = total + matrix[m][i]
        Total_row = Total_row + [total]
    return(Total_row)

def Calculate_pj_row(matrix):
    pj_row = []
    Total_row = Calculate_Total_row(matrix) 
    for x in range(0,len(Total_row)):
        pj_row = pj_row + [Total_row[x]/float(sum(Total_row))]
    return(pj_row)  
      
def Calculate_first_order(matrix):
    """"input: nij matrix (See Appendix A of JORS publication)
    output: Sum of pj"""
    pj_row = Calculate_pj_row(matrix)
    total = 0
    for m in range(0,len(pj_row)):
        total = total + pj_row[m]
    return(total)

def Calculate_Pe(matrix):
    """"input: nij matrix (See Appendix A of JORS publication)
    output: Sum of pj^2, which corresponds to pe"""
    pj_row = Calculate_pj_row(matrix)
    total = 0
    for m in range(0,len(pj_row)):
        total = total + pj_row[m]**2
    return(total)

def Calculate_third_order(matrix):
    """"input: nij matrix (See Appendix A of JORS publication)
    output: Sum of pj^3"""
    pj_row = Calculate_pj_row(matrix)
    total = 0
    for m in range(0,len(pj_row)):
        total = total + pj_row[m]**3
    return(total)
#==============================================================
#==============================================================


#=====================Fleiss kappa and variance================
#==============================================================
def Calculate_Fleiss_kappa(matrix):
    """input: nij matrix (See Appendix A of JORS publication)"""
    Pe = Calculate_Pe(matrix)
    P = Calculate_P(matrix)
    kappa = (P - Pe)/(1-Pe)
    return(kappa)

def Calculate_Fleiss_kappa_SE(n, matrix):
    """
    Input the number of raters and the nijmatrix
    Calculates variance using theory from:
    Fleiss, Joseph L., John C. Nee, and J. Richard Landis. 
    "Large sample variance of kappa in the case of different sets of raters." 
    Psychological bulletin 86.5 (1979): 974."""    
    N = len(matrix)
    #print("Number of rated subjects: " + "%i" % N)        
    Pe = Calculate_Pe(matrix)
    third_order = Calculate_third_order(matrix)
    numerator = 2*((1-Pe)**2 + 3*Pe - 1 - 2*third_order)
    denominator = N*n*(n-1)*(1-Pe)**2
    SE = (numerator/denominator)**0.5
    return(SE)
###================================================================
###================================================================


###==================Cohen kappa and variance =====================
###================================================================
def Calculate_nki(array,categories):
    nki = []
    for i in range(0,len(array[0])):
        error = 0
        nk = []
        u_total = 0
        for k in range(0,len(categories)):
            t = 0
            for l in range(0,len(array)):
              if array[l][i] == categories[k]:
                t = t + 1             #How many times did the user pick that rating? 
                u_total = u_total + 1 #Count the number of times the user made a valid rating
            nk = nk + [t]          
        for l in range(0,len(array)):
          if array[l][i] not in categories:
            error = error + 1
            u_total = u_total + 1
        if u_total != len(array):
            print("Error")
        nk = nk + [error]
        nk = [x / float(u_total) for x in nk ] #Normalize
        nki = nki + [nk]
    return(nki)

def Calculate_po(array,categories):
    po = 0
    for m in range (0, len(array)):    
      if array[m][0] == array[m][1]:
        po = po + 1
    po = po/float(len(array))
    return(po)

def Calculate_pok(array,categories):
    pok = []
    for k in range (0, len(categories)):
        t = 0
        for m in range (0, len(array)):    
            if array[m][0] == categories[k] and array[m][1] == categories[k]:
                t = t + 1
        pok = pok + [t]
    pok = [x / float(len(array)) for x in pok ] #Normalize
    return(pok)

def nki_pe(array,categories):
    nki = Calculate_nki(array,categories)
    pe = 0
    for k in range(0,len(categories)):
        pe = pe + nki[0][k]*nki[1][k]
    return(pe)        

def Calculate_Cohen_kappa(array,categories):
    po = Calculate_po(array,categories)
    pe = nki_pe(array,categories)
    kappa = (po - pe)/(1-pe)
    return(kappa)

def Calculate_Cohen_kappa_SE(array,categories):
    """Input a 2D array, ratings between two users.
    Outputs standard error of kappa. Using equation from Fleiss1979 paper"""
    matrix =  Calculate_nijmatrix(array,categories)
    SE = Calculate_Fleiss_kappa_SE(2, matrix)
    return(SE)

def Calculate_Cohen_kappa_simplisticSE(array,categories):
    """Input a 2D array, ratings between two users.
    Outputs standard error of kappa. 
    Using equation form McHugh 2012 paper. This 
    expression overestimates the standard error"""
    po = Calculate_po(array,categories)
    pe = nki_pe(array,categories)
    N = len(array)
    SE = math.sqrt(po*(1-po)/(N*(1-pe)**2))    
    return(SE)
###================================================================
###================================================================


###===================Permutated Kappa=============================
###================================================================
def Extract_two(array, r1, r2):
    """Extract two columns (raters) from a data file"""
    new = []
    N = len(array)
    for m in range(0,N):
        newrow = [array[m][r1]] + [array[m][r2]] 
        new = new + [newrow]
    return(new)    

def Keep_ValidRhatings(Dosarray, categories):
    """Input a 2D array - the ratings from two users. Only keep subjects
    where both of these raters entered correct ratings"""
    new = []
    for m in range(0,len(Dosarray)):
      if (Dosarray[m][0] in categories) and (Dosarray[m][1] in categories):  
                newrow = [Dosarray[m][0]] + [Dosarray[m][1]]
                new = new + [newrow]
    return(new)

def Calculate_permutated_kappa_tensor(array,categories):
    """Input the data file with ratings from all users.
    Function outputs the permutated-kappa tensor (PK_tensor),
    See JORS publication."""
    n = len(array[0])       #Number of raters
    PK_tensor = []
    for y in range(0,n):
      row = []
      for x in range(0,n):
            z = Extract_two(array,y,x)          #Extract two users
            z = Keep_ValidRhatings(z, categories)   #Keep only simultaneous ratings
            kappa = Calculate_Cohen_kappa(z,categories)
            SE = Calculate_Cohen_kappa_SE(z,categories)
            row = row + [[kappa, SE]]
      PK_tensor = PK_tensor + [row]
    return(PK_tensor)

def PK_tensor_per_rater(PK_tensor, rater): 
    """Input the permutated kappa tensor, and the desired rater.
    Function calculates the average between himself and the other raters.
    Outputs an array of the average and standard error"""
    n = len(PK_tensor)
    sumk = 0
    sumdiff = 0
    for x in range(0,n):                    #Sample all raters
        if x != rater:                      #Do not compare rater with herself
            sumk = sumk + PK_tensor[rater][x][0]
    aver_k = sumk/(n-1)                     #n-1, since we do not count rater with herself
    for x in range(0,n):                    #Sample all raters
        if x != rater:                      #Do not compare rater with herself
            sumdiff = sumdiff + (PK_tensor[rater][x][0]-aver_k)**2
    std_k = (sumdiff/(n-2))**0.5            #Sample standard deviation
    SE = std_k/(n-1)**0.5                  #Standard error
    return([aver_k, SE])
###================================================================
###================================================================


###=======================Output to screen=========================
###================================================================
def Fleiss_kappa_analysis(n, nij_matrix, categories):
    """Print out stats for group"""
    print("-------------------------------------------")
    print("Fleiss-kappa analysis (Group statistics):")
    print("-------------------------------------------")
    #1) Calculate probabilities for each rating
    pj = Calculate_pj_row(nij_matrix)
    #2) Calculate probability of random agreement
    Pe = Calculate_Pe(nij_matrix)
    #3) Calculate average user agreement
    P = Calculate_P(nij_matrix)
    #4) Calculate Fleiss_kappa and its standard error
    kappa = Calculate_Fleiss_kappa(nij_matrix)
    SE = Calculate_Fleiss_kappa_SE(n, nij_matrix)
    #5) Print out stats
    for k in range(0,len(categories)):
        print("Probability that '" + categories[k] + "' was chosen: " + \
              "%.3f" % pj[k])
    print("Agreement expected from chance (Pe): " + "%.3f" % Pe)
    print("Average user agreement (P): " + "%.3f" % P)
    print("kappa: " + "%.3f" % kappa)
    print("kappa Confidence interval: [" + "%.3f" % (kappa - 1.96*SE) + \
                                       ", %.3f" % (kappa + 1.96*SE) + "]\n")
    return

def Stats(array,categories, pk_tensor):
    """Print out stats for each user"""
    print("\n-------------------------------------")
    print("User statistics: ")
    N = len(array)          #How many subjects were examined?
    print("Number of rated subjects (N): " + "%i" % N)        
    n = len(array[0])       #Number of raters
    print("Number of raters (n): " + "%i" % n)        
    print("-------------------------------------")
    for user in range(0,n):
      stats = []
      u_total = 0
      print("\nUser: " + "%i" % (user))
      for j in range(0,len(categories)):      
          nxi = 0
          for i in range(0,len(array)) :
              if array[i][user] == categories[j]:
                nxi = nxi + 1         #How many times did the user 
                u_total = u_total + 1 #Count the number of times the user made a valid rating 
          stats = stats + [nxi]         
      y = user                        #Dummy variable for later use
      stats = [user / float(u_total) for user in stats ]  #Normalize by number of valid ratings
      for k in range(0,len(categories)):
          print("Probability that '" + categories[k] + "' was chosen: " + \
                "%.3f" % stats[k])
      for x in range(0,n):
          if x!=y:
            cilow = pk_tensor[x][y][0]-1.96*pk_tensor[x][y][1]
            cihigh = pk_tensor[x][y][0]+1.96*pk_tensor[x][y][1]
            print("Cohen kappa between users " + "%i" % x + " and " + \
                  "%i" % y + ": " + "%.3f" % pk_tensor[x][y][0] + \
                  " [" + "%.3f" % cilow + ", " + "%.3f" % cihigh + "]")
      pkz = PK_tensor_per_rater(pk_tensor, y)
      print("Permutated kappa average: " + "%.3f" % pkz[0])
      print("Permutated kappa standard deviation: " + "%.3f" %pkz[1])
      print("Permutated kappa confidence interval: [" \
            + "%.3f" % (pkz[0]-1.96*pkz[1]) + ", " + \
            "%.3f" % (pkz[0]+1.96*pkz[1])+"]")
    print("-------------------------------------")
    return
###================================================================
###================================================================

  
###===================Plot output graph============================
###================================================================
def Plot_PK_tensor(PK_tensor, nij_matrix, ymin, ymax, graph_filename, h, indbars):
    
    #--------------------------------------
    #Calculate Fleiss kappa
    n = len(PK_tensor)
    kappa = Calculate_Fleiss_kappa(nij_matrix)
    SE = Calculate_Fleiss_kappa_SE(n, nij_matrix)
    #--------------------------------------

    #--------------------------------------
    #Plot size and axis
    fig1 = plt.figure(figsize=(20,12))
    plt.xlim(-0.5,n-0.5)
    plt.ylim(ymin, ymax)
    ax1 = fig1.add_subplot(111)
    #--------------------------------------

    #-----------------------------------------------------------------------------------
    #Plot values from Fleiss-kappa
    plt.axhline(y=kappa + 1.96*SE, linestyle = "--", color = 'g', markersize = 30)
    plt.axhline(y=kappa - 1.96*SE, linestyle = "--", color = 'g', markersize = 30)
    plt.text(n-0.8, kappa + 1.96*SE+0.01, "$\kappa_F$ upper limit " , fontsize = 25)
    plt.text(n-0.8, kappa - 1.96*SE-0.01, "$\kappa_F$ lower limit" , fontsize = 25)
    #------------------------------------------------------------------------------------

    #------------------------------------------------------------------------
    #Build an array for the graph. There will be n(n-1)entries.
    #x-values: User #
    #y-values: permutated kappa
    #y-error: 1.96*SE
    x_values = []; y_values = []; y_error = []; counter = 0
    for i in range(0,len(PK_tensor)):
        for j in range(0,len(PK_tensor[0])):
          if i != j:
            counter = counter + 1
            x_values = x_values + [i]
            y_values = y_values + [PK_tensor[i][j][0]] 
            y_error = y_error + [1.96*PK_tensor[i][j][1]] 
    entries = len(PK_tensor)*(len(PK_tensor) - 1)
    #------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    #Plot permutated-kappas for each user
    if indbars == "no":
        ax1.errorbar(x_values, y_values, fmt='^', markersize = 12,\
                     label = "Permutated kappas")
    #If researchers wants individual error bars on each pair Cohen Kappa
    if indbars == "yes":   
        for i in range(0,entries):
            ax1.errorbar(x_values[i], y_values[i], y_error[i], \
                     elinewidth = 2, fmt = '^', markersize = 12)
    #Label each pair (i.e. [0,2])
    for i in range(0,len(PK_tensor)):
        for j in range(0,len(PK_tensor[0])):
            if i != j:
              plt.text(i+0.02, PK_tensor[i][j][0], ("(%i" % i + ", %i)" %j)\
                       , fontsize = 14)
    #--------------------------------------------------------------------------
    
    #----------------------------------------------------------
    #For each user: Plot average and CI from permutated-kappas
    x_values = []; y_values = []; y_error = []
    for r in range(0,len(PK_tensor)):
        st = PK_tensor_per_rater(PK_tensor, r)
        x_values = x_values + [r]
        y_values = y_values +  [st[0]]
        y_error = y_error + [st[1]]
    ax1.errorbar(x_values, y_values, y_error, fmt='o', elinewidth = 4, \
                 markersize = 10, color = "orange", \
                 label ="User's average of the permutated kappas")
    #-----------------------------------------------------------

    #-------------------------------------------------------------------
    #If researchers chooses to highlight a user-pair
    if h[0] == "none":
        print("\n")
    else:
        u = h.split(',')
        x = int(u[0]); y = int(u[1]);
        print("User chose to highlight the: %i," % x + "%i pair" % y)
        x_values2 = [x,y]
        y_values2 = [PK_tensor[x][y][0], PK_tensor[x][y][0]]
        y_error2 = [PK_tensor[y][x][1], PK_tensor[y][x][1]]
        ax1.errorbar(x_values2, y_values2, y_error2, fmt = 'o', \
                     elinewidth = 2, color = "black", markersize = 8)    
    #-------------------------------------------------------------------

    #-----------------------------------
    #Alert the reader of the output file
    print("Output graph saved as: %s \n" %graph_filename)    
    #-----------------------------------

    #-----------------------------------
    #Labels, legends, ticks, save graph
    #-----------------------------------
    plt.xlabel("Users", fontsize = 20)
    plt.ylabel("kappa", fontsize = 20)
    plt.legend(loc='upper left', fontsize = 20);
    plt.xticks(x_values, fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.savefig(graph_filename, format = 'jpg', dpi = 300)
###================================================================
###================================================================