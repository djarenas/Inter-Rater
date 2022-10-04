#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: October 2022
"""

import csv
import math
import matplotlib.pyplot as plt


#==============================================================
#These are the functions for reading data files
#==============================================================
def readfile_categories(filename):
    """
    The one thing it does: Categories-style files --> list
    Input: string filename. File must be tab separated.
    Output: List with k categories
    """
    rlist = []        #Lists are mutable
    with open(filename, 'r') as inF: 
      for line in csv.reader(inF, dialect="excel-tab"):
          rlist = rlist + line # "+" concatenates
    return(rlist)

def readfile_array(filename):
    """
    The one thing it does: Ratings-file --> Two-dimensional list
    Input: string filename.
    Output: A list of the form [[row1], [row2], ....]
    Filename style:    
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1

    """
    list2D = []
    with open(filename, 'r') as inF:
        for line in csv.reader(inF, dialect="excel-tab"):
            list2D = list2D + [line] 
    return(list2D)

def print2Dlist (l2):
    print("\n")
    for x in l2:
        for y in x:
            print(y, end="\t")
        print("\n")
    return


#==============================================================
#These are the functions for calculating/handling nij matrices
#==============================================================
def calc_specific_nij(i, r, list2D):
    """
    The one thing it does: Calc # of raters who chose a rating for ith item
    Inputs:
    i: item index  (int from 1 to n)
    r: rating (string )
    list2D: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1
    Output: nij (integer)
    """
    
    row = list2D[i]
    t = sum([True for j in row if j == r])
    return(t)


def calc_row_nij(i,categories, list2D):
    """
    The one thing it does: Calculates the ith row of the nij matrix 
    Inputs:
    i: item index  (int from 1 to N)
    categories: Possible ratings (list of strings)
    list2D: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1
    Output: ith row of nij matrix (list of integers)
    """
    row = []
    for c in categories:
        nij = calc_specific_nij(i,c, list2D)
        row = row + [nij]
    return(row)

def calc_nijmatrix(list2D,categories):
    """
    What it does: Calculate nij matrix. 
    nij matrix is the number of raters who chose jth category in ith item. 
    Inputs:
    categories: Possible ratings (list of strings)
    array: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1    
    Output: nij matrix ([[row1], [row2],...])
    """
    ijmatrix = []
    N = len(list2D) #Number of items that were rated
    for m in range(0,N):
        row = calc_row_nij(m,categories,list2D)
        ijmatrix = ijmatrix + [row]
    return(ijmatrix)

#==============================================================
#These functions take the nij matrix and calculate probabilities
#==============================================================
def calc_probi(ijmatrix, i):
    """
    Purpose: Calculate agreement-probability for ith item 
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: Agreement-probability for ith item (float)
    """
    k = len(ijmatrix[i])          #k: Number of categories
    n = 0.0                     #n: Number of ratings
    #How many raters evaluated ith item?
    for j in range(0,k):
        n = n + ijmatrix[i][j]
    #Calculate probability of agreement for ith item
    probi = 0.0
    for j in range(0,k):
        nij = ijmatrix[i][j]
        probi = probi + (1/(n**2-n))*(nij**2-nij)
    return(probi)

def calc_probi_column(ijmatrix):
    """
    Purpose: Calculate agreement (0 to 1) for each of the Nth items
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: Probability-agreement column (list of floats from 0 to 1) 
    """
    N = len(ijmatrix)             #N: Number of studies
    probi_column = []
    for j in range(0,N):
        probi_column = probi_column + [calc_probi(ijmatrix, j)]
    return(probi_column)

def calc_po_group(matrix):
    """
    Purpose: Obtain total agreement (0 to 1)
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: Agreement column (float from 0 to 1) 
    """
    Pi_column = calc_probi_column(matrix)
    return(sum(Pi_column)/float(len(matrix)))

def calc_total_row(matrix):
    """
    Purpose: Across all items, obtain number of total ratings for jth category 
            Not trivial since not every rater needs to evaluate every item
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: Total ratings for each category (list of integers) 
    """
    total_row = []    
    #For each item
    for j in range (0, len(matrix[0])):
        totalj = 0
        #For each category
        for i in range (0,len(matrix)):
            totalj = totalj + matrix[i][j]
        total_row = total_row + [totalj]
    return(total_row)

def calc_pj_row(matrix):
    """
    Purpose: For each jth category, obtain probability that it was chosen
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: pj_row (list of floats) 
    """
    pj_row = []
    Total_row = calc_total_row(matrix) 
    for x in range(0,len(Total_row)):
        pj_row = pj_row + [Total_row[x]/float(sum(Total_row))]
    return(pj_row)  
      
def calc_first_order(matrix):
    """
    Purpose: Sum of first orders of pj (Should always equal one)
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: float (0 to 1)
    """    
    pj_row = calc_pj_row(matrix)
    return(sum(pj_row))

def calc_pe(matrix):
    """
    Purpose: Calculate agreement by chance: Sum of second orders of pj
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: float (0 to 1)
    """    
    pj_row = calc_pj_row(matrix)
    pe = sum([x**2 for x in pj_row])
    # if (pe == 1):
    #         print("Warning: A probability of agreement by chance \
    #               was calculated to be one. This means only one rating was \
    #               chosen by n > 1 users in all items.")
    return(pe)

def calc_third_order(matrix):
    """
    Purpose: Sum of third orders of pj
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: float (0 to 1)
    """    
    pj_row = calc_pj_row(matrix)
    return(sum([x**3 for x in pj_row]))


#================================================================
#These functions calculate Fleiss kappa and its variance (n > 2)
#================================================================
def calc_fleiss_kappa(matrix):
    """
    Purpose: Take the nij matrix and calculate Fleiss kappa (groups > 2)
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: float (0 to 1)
    """    
    Pe = calc_pe(matrix)
    P = calc_po_group(matrix)
    kappa = (P - Pe)/(1-Pe)
    return(kappa)

def calc_fleiss_kappa_SE(n, matrix):
    """
    Purpose: Calculate variance of Fleiss kappa (groups > 2)
             Based on an equation from Fleiss 1979 manuscript
    Inputs:
    matrix: nij matrix (list in the format [[row1], [row2],...], )
    Output: Square root of variance (float)
    """ 
    N = len(matrix)
    first_order = calc_first_order(matrix)
    if abs(first_order - 1) > 0.01:
        print ("Your probabilities are not adding to one")
    Pe = calc_pe(matrix)
    third_order = calc_third_order(matrix)
    numerator = 2*((1-Pe)**2 + 3*Pe - 2*third_order - 1)
    denominator = N*n*(n-1)*(1-Pe)**2
    SE = (numerator/denominator)**0.5
    return(SE)


#================================================================
#These functions calculate Cohen kappa and its variance (n = 2)
#================================================================
def calc_cohen_nij(array,categories):
    """
    Purpose: Return nij matrix for n = 2 (All items must have valid ratings)
    Inputs:
    categories: Possible ratings (list of strings)
    array: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1    
    Output: nij matrix ([[row1], [row2],...])
    """
    #---------------------------------------------
    #Check that all items have two valid ratings
    #---------------------------------------------
    error = 0
    for i in range(0,len(array[0])):
        for l in range(0,len(array)):
            if array[l][i] not in categories:
                error = error + 1
    if (error > 0):
        print("There were ", error, "items without two valid ratings")
    nijmatrix = calc_nijmatrix(array,categories)
    return(nijmatrix)

def calc_cohen_kappa(array,categories):
    """
    Purpose: Return Cohen kappa for n = 2 (All items must have valid ratings)
    Inputs:
    categories: Possible ratings (list of strings)
    array: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1    
    Output: nij matrix ([[row1], [row2],...])
    """
    matrix = calc_cohen_nij(array,categories)
    pe = calc_pe(matrix)
    po = calc_po_group(matrix)
    kappa = (po - pe)/(1-pe)
    return(kappa)

def calc_cohen_kappa_SE(array,categories):
    """
    Purpose: Calculates variance for Cohen kappa (n = 2) 
    Inputs:
        categories: Possible ratings (list of strings)
        array: List in the format [[row1], [row2],...], where 
                     rater 1   rater 2   rater 3
            study 1    0         0          1
            study 2    -1        1          1    
    Output: Standard deviation of Cohen kappa (float)
    """
    matrix =  calc_cohen_nij(array,categories)
    SE = calc_fleiss_kappa_SE(2, matrix)
    return(SE)

def calc_cohen_kappa_simplisticSE(array,categories):
    """
    Purpose: Calculates variance for Cohen kappa (n = 2) 
    Using equation form McHugh 2012 paper. 
    This expression overestimates the standard error
    
    Inputs:
        categories: Possible ratings (list of strings)
        array: List in the format [[row1], [row2],...], where 
                     rater 1   rater 2   rater 3
            study 1    0         0          1
            study 2    -1        1          1    
    Output: Standard deviation of Cohen kappa (float)
    """
    matrix = calc_cohen_nij(array,categories)
    pe = calc_pe(matrix)
    po = calc_po_group(matrix)
    N = len(array)
    SE = math.sqrt(po*(1-po)/(N*(1-pe)**2))    
    return(SE)

#================================================================
#These functions permutate pairs of raters and find the kappas
#================================================================
def extract_two(list2D, r1, r2):
    """
    Purpose: Extracts two raters (indices r1 and r2) from an n > 2 group
    Inputs:
        array: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
            study 1    0         0          1
            study 2    -1        1          1    
    Output: Array of two users as a 2D list ([r11, r12], [r21, r22], ...)
    """    
    new = []
    N = len(list2D)
    for j in range(0,N):
        newrow = [list2D[j][r1]] + [list2D[j][r2]] 
        new = new + [newrow]
    return(new)    

def keep_ValidRatings(list2D, categories):
    """
    Purpose: In array for two raters, keep only items with two valid ratings
    Inputs:
        list2D: ([r11, r12], [r21, r22], ...)
    Output: Array of two users in the format ([r11, r12], [r21, r22], ...)
    """        
    new = []
    N = len(list2D)
    for i in range(0,N):
      if (list2D[i][0] in categories) and (list2D[i][1] in categories):  
                newrow = [list2D[i][0]] + [list2D[i][1]]
                new = new + [newrow]
    return(new)

def calc_permutated_kappa_tensor(array,categories):
    """
    Purpose: Calculates permutated kappa tensor    
    Inputs:
        categories: Possible ratings (list of strings)
        array: List in the format [[row1], [row2],...], where 
                     rater 1   rater 2   rater 3
            study 1    0         0          1
            study 2    -1        1          1    
    Output: 
    Returns a tensor where for each user the kappa and the standard deviation 
    are given: PK_tensor[raterx][ratery][kappa,SE]
    List in the form: [ [ [kappa12, SE12], [kappa13, SE13], ...], ... )
    """    
    n = len(array[0])       #Number of raters
    pk_tensor = []
    #print("Calculating permutated kappa tensor")
    for y in range(0,n):
        row = []
        for x in range(0,n):
            z = extract_two(array,y,x)          #Extract two users
            z = keep_ValidRatings(z, categories)   #Keep only simultaneous ratings
            #Comparing the rater to herself is trivial
            if (x == y):
                kappa = 1
                SE = 0
            else:            
                #Make sure that pe is not equal to one
                ijmatrix = calc_nijmatrix(z, categories)
                pe = calc_pe(ijmatrix) #Check that pe is not equal to one
                if (pe == 1): 
                    print("Error in users: ", str(x), ", ", str(y) )
                #Calculate kappa, SE if pe is not equal to one
                kappa = calc_cohen_kappa(z,categories)
                SE = calc_cohen_kappa_SE(z,categories)
            #Update the row
            row = row + [[kappa, SE]]
        #Update the tensor 
        pk_tensor = pk_tensor + [row]
    return(pk_tensor)

def calc_pk_tensor_average_rater(pk_tensor, rater): 
    """
    Purpose: Average cohen kappa for rater across pair permutations    
    Input: 
        PK_tensor[raterx][ratery][kappa,SE] (multidimentional list)
        rater: Which rater you want to analyze (integer)
    Output: Cohen kappa's [Average, standard deviation]
    """
    n = len(pk_tensor)
    sumk = 0
    sumdiff = 0
    for x in range(0,n):                    #Sample all raters
        if x != rater:                      #Do not compare rater with herself
            sumk = sumk + pk_tensor[rater][x][0]
    aver_k = sumk/(n-1)                     #n-1, since we do not count rater with herself
    for x in range(0,n):                    #Sample all raters
        if x != rater:                      #Do not compare rater with herself
            sumdiff = sumdiff + (pk_tensor[rater][x][0]-aver_k)**2
    std_k = (sumdiff/(n-2))**0.5            #Sample standard deviation
    SE = std_k/(n-1)**0.5                  #Standard error
    return([aver_k, SE])


#================================================================
#These functions are for printing and plotting the results
#================================================================
def fleiss_kappa_report(n, nij_matrix, categories):
    """
    Purpose: This function calculates the user agreement statistics (fleiss)
    and outputs them to screen.
    Input: 
    n: # of raters (integer)
    nij_matrix: (list)
    categories: (list of strings)
    """
    pj = calc_pj_row(nij_matrix)   
    Pe = calc_pe(nij_matrix)       #Prob of agreement by chance
    P = calc_po_group(nij_matrix)  #Agreement (0 to 1)
    kappa = calc_fleiss_kappa(nij_matrix)
    SE = calc_fleiss_kappa_SE(n, nij_matrix)

    print("-------------------------------------------")
    print("Fleiss-kappa analysis (Group statistics):")
    print("-------------------------------------------")
    for k in range(0,len(categories)):
        print("Probability that '" + categories[k] + "' was chosen: " + \
              "%.3f" % pj[k])
    print("Agreement expected from chance (Pe): " + "%.3f" % Pe)
    print("Average user agreement (P): " + "%.3f" % P)
    print("kappa: " + "%.3f" % kappa)
    print("kappa Confidence interval: [" + "%.3f" % (kappa - 1.96*SE) + \
                                       ", %.3f" % (kappa + 1.96*SE) + "]\n")
    return

def rater_report(array,categories, pk_tensor):
    """
    Purpose: This function reports the results for each rater 
    Input: 
    categories: Possible ratings (list of strings)
    array: List in the format [[row1], [row2],...], where 
                 rater 1   rater 2   rater 3
        study 1    0         0          1
        study 2    -1        1          1   
    pk_tensor: [ [ [kappa12, SE12], [kappa13, SE13], ...], ... )
    """
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
      #Count the number of absent or invalid ratings
      none = 0 
      for i in range(0,len(array)):
          if array[i][user] not in categories:
              none = none + 1
      y = user                        #Dummy variable for later use
      stats = [user / float(u_total) for user in stats ]  #Normalize by number of valid ratings

      #Print out the stats
      if none > 0:
            print("Number of absent or invalid ratings: " + \
                "%i" % none + ". Proportion of total: %.2f" % float(none/N))

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
      pkz = calc_pk_tensor_average_rater(pk_tensor, y)
      print("Permutated kappa average: " + "%.3f" % pkz[0])
      print("Permutated kappa standard deviation: " + "%.3f" %pkz[1])
      print("Permutated kappa confidence interval: [" \
            + "%.3f" % (pkz[0]-1.96*pkz[1]) + ", " + \
            "%.3f" % (pkz[0]+1.96*pkz[1])+"]")
    print("-------------------------------------")
    return

def plot_PK_tensor(pk_tensor, nij_matrix, ymin, 
                   ymax, graph_filename, h, indbars):
    """
    Purpose: Saves the results into a graph (jpg file)
    Input: 
    pk_tensor: [ [ [kappa12, SE12], [kappa13, SE13], ...], ... )
    nij_matrix: (list)
    ymin, ymax: Minimum and maximum of y axes (float)
    graph_filename: (string)
    highlight: Pairs you wish to highlight (list of integers)
    indbars: Option for individual error bars on each cohen kappa ("yes or no")
    """    
    #--------------------------------------
    #Calculate Fleiss kappa
    n = len(pk_tensor)  #Number of raters
    kappa = calc_fleiss_kappa(nij_matrix)
    SE = calc_fleiss_kappa_SE(n, nij_matrix)
    #--------------------------------------

    #--------------------------------------
    #Plot size and axis
    fig1 = plt.figure(figsize=(20,12))
    plt.xlim(-1,n-0.5)
    plt.ylim(ymin, ymax)
    ax1 = fig1.add_subplot(111)
    #--------------------------------------

    #-----------------------------------------------------------------------------------
    #Plot values from Fleiss-kappa
    plt.axhline(y=kappa + 1.96*SE, linestyle = "--", color = 'g', markersize = 20)
    plt.axhline(y=kappa - 1.96*SE, linestyle = "--", color = 'g', markersize = 20)
    plt.text(-0.95, kappa + 2.2*SE, "$\kappa_F$ upper limit " , fontsize = 25)
    plt.text(-0.95, kappa - 2.2*SE-0.02, "$\kappa_F$ lower limit" , fontsize = 25)
    #------------------------------------------------------------------------------------

    #------------------------------------------------------------------------
    #Build an array for the graph. There will be n(n-1)entries.
    #x-values: User #
    #y-values: permutated kappa
    #y-error: 1.96*SE
    x_values = []; y_values = []; y_error = []; counter = 0
    for i in range(0,len(pk_tensor)):
        for j in range(0,len(pk_tensor[0])):
          if i != j:
            counter = counter + 1
            x_values = x_values + [i]
            y_values = y_values + [pk_tensor[i][j][0]] 
            y_error = y_error + [1.96*pk_tensor[i][j][1]] 
    entries = len(pk_tensor)*(len(pk_tensor) - 1)
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
    for i in range(0,len(pk_tensor)):
        for j in range(0,len(pk_tensor[0])):
            if j < len(pk_tensor[0]) - 1:  
              delta = pk_tensor[i][j][0] - pk_tensor[i][j+1][0]
            if i != j:
              plt.text(i+0.02, pk_tensor[i][j][0], ("(%i" % i + ", %i)" %j)\
                       , fontsize = 14)
    #--------------------------------------------------------------------------
    
    #----------------------------------------------------------
    #For each user: Plot average and CI from permutated-kappas
    x_values = []; y_values = []; y_error = []
    for r in range(0,len(pk_tensor)):
        st = calc_pk_tensor_average_rater(pk_tensor, r)
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
        y_values2 = [pk_tensor[x][y][0], pk_tensor[x][y][0]]
        y_error2 = [pk_tensor[y][x][1], pk_tensor[y][x][1]]
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
    plt.legend(loc='upper right', fontsize = 20);
    plt.xticks(x_values, fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.savefig(graph_filename, format = 'jpg', dpi = 200)
###================================================================
###================================================================