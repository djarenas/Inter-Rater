#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:35:52 2017
@author: daniel
Last update: January 2023
"""

import numpy as np
import matplotlib.pyplot as plt

#==============================================================
#Functions for reading the data files
#==============================================================
def read_file(filename):
    """
    Purpose: Tab-separated file -> Numpy-array of strings
    
    Parameters:
    ----------
    filename : string.
        Include the extension (i.e. '.txt') if the file has one.
        The file you want to read must be tab-separated, and have
        a equal number of columns in every row.

    Returns: 
    -------
        Numpy array of strings.
    
    """
    
    ratings_array = np.loadtxt(filename, dtype = str, delimiter = '\t')
    
    return(ratings_array)


def check_input_data(categories, ratings_array):
    """
    Purpose: Check that the input data can be used in an inter-rater analysis.

    Parameters
    ----------
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.
        
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user
                         rater 1   rater 2   rater 3
        item 1    0         0          1
        item 2    -1        1          1

    Returns
    -------
    bool
        True if the data is suitable for inter-rater analysis.
        False if errors were found that would make the data unusable for 
        inter-rater analysis.
    """
    
    #Check that there are at least two categories.
    #Otherwise: Throw error and terminate.
    if len(categories) < 2:
        print("Error in categories file: \
              There should at least be 2 valid categories. Terminating.")
        return False
    
    #Check that there are at least two evaluators in the ratings data.
    #Otherwise: Throw error and terminate.
    if len(ratings_array[0]) < 2:
        print("Error in data file: \
              There should data from at least 2 raters. Terminating.")
        return False
    
    #Check that there are at least 10 items evaluated.
    #Otherwise, throw a warning but do not terminate.
    if (len(ratings_array) < 12):
        print("Warning for data file. Less than 10 items were evaluated. \
              Will continue; but, be aware that this number may be too low for \
              the equations used in the estimation of the variances")

    return True


#==============================================================
#Functions for calculating/handling nij matrices
#==============================================================


def calc_specific_nij(i, k, ratings_array):
    """
    Purpose: Calculate # of raters who chose kth category for ith item

    Parameters
    ----------
    i : Integer
        Index for evaluated item
    r : String
        Specific rating (ratings)
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user,
        Example:
                         rater 1   rater 2   rater 3
        item 1    0         0          1
        item 2    -1        1          1
    
    Returns
    -------
        Integer
    """
    
    return (sum(ratings_array[i] == k))


def calc_nijmatrix(ratings_array,categories):
    """
    Purpose: Calculate the nij 2D array. 
    The output will keep the original order of categories specified in
    the original categories file.

    Parameters
    ----------
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user
                         rater 1   rater 2   rater 3
        item 1    0         0          1
        item 2    -1        1          1
        
        This ratings array must contain at least two columns.
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.
        Must contain at least two elements.

    Returns 
    -------
    2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item. 
    """    
    
    #Calculate the first two columns. 
    ni1 = np.sum(ratings_array == categories[0], axis = 1)
    ni2 = np.sum(ratings_array == categories[1], axis = 1)
    nij = np.dstack((ni1,ni2))

    #Continue for the rest of the categories and append them to the nij matrix.
    for k in range(2,len(categories)):
        nik = np.sum(ratings_array == categories[k], axis = 1)
        nij = np.dstack((nij, nik))
    
    #Return a numpy array with nth rows (number of evaluated items) 
    #and kth columns (number of categories)
    return(nij[0])


def calc_Probi_column(nij_array):
    """
    Purpose: Calculate agreement (0 to 1) for each of the N items.
    Equation: Pi = 1/(n^2-1)*sumj(nij^2 - nij)
    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.
    Returns 
    -------
    1D Numpy array of floats
        Probability-agreement for each item
    """
    
    #Number of evaluators per item.
    n = np.sum(nij_array, axis = 1)
    
    #Use of Lambda expression to calculate sum(nij^2 - nij)
    agr = lambda x: x*(x-1)
    prob_i = np.sum(agr(nij_array),axis =1)
    
    #Normalize over n(n-1)
    prob_i = (prob_i)/(n*(n-1))
    
    #Return 1D array of float
    return prob_i


def calc_Po(nij_array):
    """
    Purpose: Calculate the total agreement (0 to 1)

    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    Float
        Total agreement probability (0 to 1).

    """
    prob_i = calc_Probi_column(nij_array)
    return (sum(prob_i)/len(prob_i))


def calc_Pj_row(nij_array):
    """
    Purpose: Calculate the probability each jth category was chosen.

    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    1D array of floats
        The probability that each jth category was chosen.
    """
    
    #Total number of evaluations
    #Not trivial since not every rater has to evaluate every item.
    n = sum(np.sum(nij_array, axis = 1))
    
    return sum(nij_array)/n


def calc_Pe(nij_array):
    """
    Purpose: Calculate the probability that the raters agree by random chance.
    This is the second moment of the Pj row.

    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    pe : Float (0 to 1)
        Probability that raters agree by random chance.

    """
    pj_row = calc_Pj_row(nij_array)
    
    #The random chance of agreement is the second order of Pj 
    pe = np.sum(pj_row**2)
        
    return pe
    

def check_nijmatrix(nij_array):
    """
    Purpose: Check that each item had enough valid ratings to 
    perform inter-rater analysis.
    
    Also check if raters chose one category for every single item. 
    This leads to random agreement chance of 1, and trivial inter-rater analysis.

    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    Boolean
        False if inter-rater analysis cannot be performed.
    """
    
    #Number of evaluators per item
    n = np.sum(nij_array, axis = 1)

    #Check each item had at least 2 evaluators.
    #If so, return false and end the function.
    if (sum(n < 2) != 0):
        print("Error. Some items did not have at least 2 valid ratings")
        itemindex = np.where(n < 2)
        print("Array indices:")
        print(itemindex)
        return False
    
    #Check that random chance of agreement is not equal to one.
    pe = calc_Pe(nij_array)
    if (pe == 1):
        print("Random agreement was calculated to be 1")
        return False

    #Otherwise...
    return True


def calc_Kappa(nij_array):
    """
    Purpose: Take the nij 2D array and calculate the kappa
    
    If (n > 2): "Group". Yields Fleiss kappa. 
    Designed to work even for situations
    where not every rater rated every item. However, each item must have
    been rated by at least two raters.
    
    If (n = 2) Function yields Cohen kappa if the array has only two raters
    and ALL the items were rated by the two raters.
    
    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    Float (-infinity to 1)
        Kappa 
    """

    #Calculate the average Fleiss Kappa
    po = calc_Po(nij_array)
    pe = calc_Pe(nij_array) 
    kappa = (po - pe)/(1-pe)
    
    return kappa


def estimate_KappaVariance(nij_array):
    """
    Purpose: Calculate the variance of the kappa statistic.
    The equation for variance comes from Fleiss 1979 manuscript:
    "Large Sample Variance of Kappa in the Case of Different Sets of Raters".
    Equation [12] was rearranged so that we only have to calculate Pe
    and the third order of the Pj row. 
    One modification/approximation is that since not every rater 
    evaluates each item, the numbers of raters as approximated 
    as the average of ratings per item.
    
    Parameters
    ----------
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.

    Returns
    -------
    Float
        Variance of kappa
    """
    
    #Number of items
    N = len(nij_array)
    
    #Number of raters (approximated to account for the fact that not
    #every rater rates every item)
    n = np.average(np.sum(nij_array, axis = 1))   
    
    #Orders of Pj row
    pj = calc_Pj_row(nij_array)
    pe = calc_Pe(nij_array)
    third_order = sum(pj**3)

    #Rearranged Equation 12 in Fleiss1979 manuscript
    numerator = 2*((1-pe)**2 + 3*pe - 2*third_order - 1)
    denominator = N*n*(n-1)*(1-pe)**2
    
    return (numerator/denominator)


def extract_TwoRaters(ratings_array, categories, r1, r2):
    """
    Purpose: 
    Return a subarray containing only data from two specific raters and
    only with items for which both users gave valid ratings.

    Parameters
    ----------
    ratings_array : 2D Numpy array of strings.
    r1 : integer
        First rater you want to select.
    r2 : integer
        Second rater you want to select

    Returns
    -------
    2D numpy array
    """
    #Shallow copy
    sc = ratings_array[:, [r1,r2]]
    
    #Deep copy
    dc = sc.copy()
    
    #Keep only rows for which each element can be found in categories
    r = dc[np.all(np.isin(dc, categories), axis = 1)]
    
    return r


def calc_PermutatedKappa(ratings_array, categories):
    """
    Purpose: Calculate the kappa for each permutation of rater pairing. Return
    a 2D array where [x][y] yields the kappa beetween raters x and y.

    Parameters
    ----------
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user
                         rater 1   rater 2   rater 3
        item 1    0         0          1
        item 2    -1        1          1
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.

    Returns
    -------
        Tuple of two numpy parrays
            First ratings_array: 2D array of permutated kappas
            Second ratings_array: 2D array of the variance of the kappas.

    """

    #Number of raters
    n = len(ratings_array[0])
    
    #Initialize the 2D arrays for the kappa and the variance with zeros.
    permutatedkappa = np.zeros((n,n))
    permutatedkvariance = np.zeros((n,n))
    
    #For every rater
    for x in range(0,n):
        #Compared to every other rater
        for y in range(0,n):
            #Except themselves
            if (x != y):
                
                #Extract data on this raters-pair, and only the items
                #where both raters gave valid ratings.
                subarray = extract_TwoRaters(ratings_array, categories, x, y)
            
                nijxy = calc_nijmatrix(subarray,categories)
                
                #If the nij matrix  is usable...
                if (check_nijmatrix(nijxy)):
                    permutatedkappa[x][y] = calc_Kappa(nijxy)
                    permutatedkvariance[x][y] = estimate_KappaVariance(nijxy)
                
                #if the matrix is not usable...
                else:
                    print("Error in calculating permutated kappa tensor")
                    print("The Error occured when comparing these two users:")
                    print(x, y)
            
            #Kappa of rater with self is trivial
            else:
                permutatedkappa[x][y] = 1
                permutatedkvariance[x][y] = 0
        #Done comparing with every other rater
    
    #Done with every rater
    
    return (permutatedkappa, permutatedkvariance)


def average_PermutatedKappaPerRater(permutatedkappa, r):
    """

    Parameters
    ----------
    permutatedkappa: Tuple of two numpy array of floats.
            First ratings_array: 2D array of permutated kappas
            Second ratings_array: 2D array of the variance of the kappas.
    r : integer
        Specific rater

    Returns
    -------
    Tuple of two integers
        (Average of permutated kappa for this user, variance)

    """
    
    n = len(permutatedkappa)
    
    #Get rater row
    b = permutatedkappa[0][r]
    b = np.delete(b,r)


    #Average, std, standard error    
    #Standard deviation calculated by using unbiased estimator
    av = np.average(b)
    std = np.std(b,ddof=1)
    se = std/len(b)**0.5
    
    return(av, se)


def print_GroupReport(n, nij_array, categories):
    """
    Purpose: Print the ratings and group kappa statistics.

    Parameters
    ----------
    n : Integer
        Number of raters
    nij_array: 2D Numpy array of integers.
        Each element equal to # of raters who chose jth category in ith item.
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.
        Must contain at least two elements.

    Returns
    -------
    None.

    """
    
    #Calculate values
    pj = calc_Pj_row(nij_array)
    pe = calc_Pe(nij_array)       
    po = calc_Po(nij_array)
    fkappa = calc_Kappa(nij_array)
    fkv = estimate_KappaVariance(nij_array)

    print("-------------------------------------------")
    print("Fleiss-kappa analysis (Group statistics):")
    print("-------------------------------------------")
    for k in range(0,len(categories)):
        print("Probability that '" + categories[k] + "' was chosen: " \
              + "%.3f" % pj[k])
    print("Agreement expected from chance (Pe): " + "%.3f" % pe)
    print("Average user agreement (P): " + "%.3f" % po)
    print("kappa: " + "%.3f" % fkappa)
    print("kappa Confidence interval: [" + "%.3f" % (fkappa - 1.96*fkv**0.5) + \
                                       ", %.3f" % (fkappa + 1.96*fkv**0.5) + "]\n")
    
    return


def get_RatingsHistogram (ratings_array, categories, r):
    """
    Purpose: For a specific rater, calculate the number of valid ratings
    and how many times they assign a specific rating. 

    Parameters
    ----------
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.
    r : Integer
        Specific rater

    Returns
    -------
    results : List
        [Total counts, valid counts, list of counts]
        where "list of counts" has the same order as in the input "categories"

    """
    
    #Get column specific to this rater
    rater_specific = ratings_array[:,r]
    
    #Count 
    total = len(rater_specific)
    
    #Count the valid ratings
    valids_boolean = np.isin(rater_specific, categories)
    valids = sum(valids_boolean)
    
    #Keeps same order as the one specified in categories
    counts = []
    for k in categories:
        counts.append(sum(rater_specific == k))
        
    #Tuple
    results = (total, valids, counts)
    
    return results


def print_RaterReports(ratings_array,categories, permukappa):
    """
    Purpose: Print the ratings and kappa statistics for each rater.

    Parameters
    ----------
    ratings_array : 2D Numpy array of strings.
        Each row contains one item, and each column a rating from each user
    categories : 1D Numpy array of strings.
        Each string corresponds to one of the valid ratings for inter-rater
        analysis.    
    permukappa: Tuple of two numpy array of floats.
            First ratings_array: 2D array of permutated kappas
            Second ratings_array: 2D array of the variance of the kappas.
    R
    eturns
    -------
    None.

    """
    #Calculate values
    N = len(ratings_array)          #How many items were examined?
    n = len(ratings_array[0])       #Number of raters  

    
    print("\n-------------------------------------")
    print("User statistics: ")
    print("Number of rated subjects (N): " + "%i" % N)        
    print("Number of raters (n): " + "%i" % n)        
    print("-------------------------------------")
    
    #For every rater
    for rater in range(0,n):

        print("\nRater: " + "%i" % (rater))
      
        #Calculate number of valid and invalid ratings
        ratings = get_RatingsHistogram(ratings_array, categories, rater)
        invalids = ratings[0]-ratings[1]
              
        #Print out the stats
      
        #Invalid ratings
        if invalids > 0:
            print("Number of absent or invalid ratings: " + \
                  "%i" % invalids + \
                  ". Proportion of total: %.2f" % float(invalids/ratings[0]))

        #Print out every time each category was chosen
        for i in range(0,len(categories)):
            print("Probability that '" + categories[i] + "' was chosen: " + \
                 "%.3f" % float(ratings[2][i]/ratings[1]))

        #Print out the kappas between the rater pairs
        for y in range(0,n):
            x = rater
            
            #Do not print out agreement with self
            if x!=y:
                cilow = permukappa[0][x][y] - 1.96*permukappa[1][x][y]**0.5
                cihigh = permukappa[0][x][y] + 1.96*permukappa[1][x][y]**0.5
                print("Cohen kappa between users " + "%i" % y + " and " + \
                  "%i" % x + ": " + "%.3f" % permukappa[0][x][y] + \
                  " [" + "%.3f" % cilow + ", " + "%.3f" % cihigh + "]")
        
        #Print average of permutated kappas
        pka = average_PermutatedKappaPerRater(permukappa, rater)
    
        print("Permutated kappa average: " + "%.3f" % pka[0])
        print("Permutated kappa standard Error: " + "%.3f" %pka[1])
        print("Permutated kappa confidence interval: [" \
              + "%.3f" % (pka[0]-1.96*pka[1]) + ", " + \
                  "%.3f" % (pka[0]+1.96*pka[1])+"]")
        print("--------------")
      
        
def plot_Kappas(fleisskappa, permutatedkappa, \
                graph_filename = "Figure.jpg", ylimits = (0,1)):
    """
    Purpose: Plot and save a jpg file visualizing the group and permutated
    kappas.

    Parameters
    ----------
    fleisskappa : Tuple of two floats.
        (Fleiss kappa, Fleiss kappa variance)
    permutatedkappa : Tuple of two 2D Numpy Arrays
            First ratings_array: 2D array of permutated kappas
            Second ratings_array: 2D array of the variance of the kappas.
    graph_filename : string, optional
        File name for the figure file. The default is "Figure.jpg".
    ylimits : Tuple of two floats, optional
        (ymin, ymax). The default is (0,1).

    Returns
    -------
    None.

    """
    #-------------------------------------------------------------------------
    #Prepare figure size and axes
    #-------------------------------------------------------------------------

    #Number of raters
    n = len(permutatedkappa[0])  

    #Figure size and axes    
    figsize = (20,12)
    plt.figure(figsize=figsize)
    plt.xlim(-0.5,n-0.5)
    plt.ylim(ylimits[0], ylimits[1])
    #ax1 = fig1.add_subplot(111)

    #-------------------------------------------------------------------------        
    #Plot Fleiss kappa confidence intervals
    #-------------------------------------------------------------------------

    kappa = fleisskappa[0]
    sd = fleisskappa[1]
    #Plot horizontal lines at the upper and lower limit of confidence interval
    plt.axhline(y=kappa + 1.96*sd, linestyle = "--", c = 'g', markersize = 20)
    plt.axhline(y=kappa - 1.96*sd, linestyle = "--", c = 'g', markersize = 20)    
    #Optional: Put text labeling these are the upper and lower limits
    #plt.text(-0.95, kappa + 2.2*sd, "$\kappa_F$ upper limit " , fontsize = 25)
    #plt.text(-0.95, kappa - 2.2*sd, "$\kappa_F$ lower limit" , fontsize = 25)

    #-------------------------------------------------------------------------
    #Plot kappa for each permutation of pairs
    #-------------------------------------------------------------------------

    #From input parameter: First array [0] is the averages of the kappas
    pkappa = permutatedkappa[0]    
    
    #For each rater...
    for x in range(0,n):
        #Compared to other raters...
        for y in range(0,n):
            #except itself...
            if x != y:
                #Plot the value for x,y kappa
                plt.errorbar(x, pkappa[x][y], fmt='^', c = 'blue', markersize = 12,\
                     label = "")
                #Plot text labeling the x,y rater pairs
                plt.text(x+1/n/4, pkappa[x][y], f"({x},{y})", fontsize = 14)

    #-------------------------------------------------------------------------
    #Plot the average of the permutated kappas
    #-------------------------------------------------------------------------

    #For each rater
    for x in range(0,n):
        #Calculate the statistics of the permutated kappa array
        av_pkappa = average_PermutatedKappaPerRater(permutatedkappa, x)
        
        #Plot the average and confidence interval
        plt.errorbar(x, av_pkappa[0], av_pkappa[1], fmt='o', elinewidth = 4, \
                      markersize = 10, color = "orange", \
                      label ="")

    #-------------------------------------------------------------------------
    #Save the plot to file
    #-------------------------------------------------------------------------
    
    plt.xlabel("Raters", fontsize = 20)
    plt.ylabel("kappa", fontsize = 20)
    plt.legend(loc='upper right', fontsize = 20);
    plt.xticks(range(0,n), fontsize = 20)
    #plt.yticks(fontsize = 20)
    plt.savefig(graph_filename, format = 'jpg', dpi = 200)
    # #Alert the reader of the output file
    # print("Output graph saved as: %s \n" %graph_filename)    
    # #-----------------------------------

    #-----------------------------------
    #Labels, legends, ticks, save graph
    #-----------------------------------
