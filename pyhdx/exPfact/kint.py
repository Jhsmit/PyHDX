#!/usr/bin/python

# This script calculates the intrinsic deuteration rates of each residue 
# of a protein into D2O
# The argument has to be the sequence of the protein in one letter code
# The temperature (in Kelvin) and the pH can be changed
# This script is based on the work of Bai and al. (1993) and the excel table proposed on the Englander group website. Notice that the parameters for E and D has been updated to take into account the work of Mori and al. (1997)
# kint is in min-1


from math import *
from pyhdx.exPfact.constants import *
import numpy as np



def calculate_kint_for_sequence(first_residue, last_residue, seq, temperature, pH):
    prolines = []
    kint = np.zeros((last_residue))
    kint.fill(-1)
    res1=""
    jj = 0
    for assignment in range(1, len(seq)+1):
        res=seq[jj]
        if not res1 == "":
           if assignment - first_residue == 0:
               kint[assignment-1] = -1
           elif seq[jj] == "P":
               kint[assignment-1] = -1
               prolines.append(first_residue + jj)
           else:
               kint[assignment-1] = calculate_kint_per_residue(res1,res,assignment,len(seq), temperature, pH)
        print("**",assignment,len(seq))
        print("***",kint[assignment-1])
        jj += 1
        res1=res
    print("Residue\tkint")
    for residue, value in zip([x for x in range(1, last_residue+1)], kint):
        print("{}\t{}".format(residue, value))

    return kint, prolines

# This function calculate the kint of a residue
# The first argument is the residue i and the second the residue i-1

def calculate_kint_per_residue(residue1, residue2, num, length, temperature, pH):

    lamb1 = acid(residue2, temperature, pH, "lamb")
    rho1 = acid(residue1, temperature, pH, "rho")
    
    if num == 2: rho1 += rho_Nterm_acid
    elif num == length: lamb1=lamb1+lamb_Cterm_acid
    
    Fa = 10**(lamb1+rho1)
    lamb2 = base(residue2, temperature, pH, "lamb")
    rho2 = base(residue1, temperature, pH, "rho")
    
    if num == 2: rho2 += +rho_Nterm_base
    elif  num == length: lamb2=lamb2+lamb_Cterm_base

    Fb = 10**(lamb2+rho2)
    
    kint = Fa * ka * get_D(pH) * get_Fta(temperature) * 60 + Fb * kb * get_OD(pH) * get_Ftb(temperature) * 60 + Fb * \
           kw * get_Ftw(temperature) * 60

    return kint	


def acid(residue, temperature, pH, value):

    if residue == "H":
        lamb = log10(10**(-0.8-pH)/(10**(-get_pK_his(temperature))+10**(-pH))+10**(-get_pK_his(temperature))/(10**(-get_pK_his(temperature))+10**(-pH)))
        rho = log10(10**(-0.51-pH)/(10**(-get_pK_his(temperature))+10**(-pH))+10**(-get_pK_his(temperature))/(10**(-get_pK_his(temperature))+10**(-pH)))			
    elif residue == "D":
        lamb = log10(10**(-0.9-pH)/(10**(-get_pK_asp(temperature))+10**(-pH))+10**(0.9-get_pK_asp(temperature))/(10**(-get_pK_asp(temperature))+10**(-pH)))			    
        rho = log10(10**(-0.12-pH)/(10**(-get_pK_asp(temperature))+10**(-pH))+10**(0.58-get_pK_asp(temperature))/(10**(-get_pK_asp(temperature))+10**(-pH)))			   	
    elif residue == "E":
        lamb = log10(10**(-0.6-pH)/(10**(-get_pK_glu(temperature))+10**(-pH))+10**(-0.9-get_pK_glu(temperature))/(10**(-get_pK_glu(temperature))+10**(-pH)))			    
        rho = log10(10**(-0.27-pH)/(10**(-get_pK_glu(temperature))+10**(-pH))+10**(0.31-get_pK_glu(temperature))/(10**(-get_pK_glu(temperature))+10**(-pH)))			   	
    else:
        lamb = para[residue][0]
        rho = para[residue][1]
    if value == "lamb":
        return lamb
    elif value == "rho":
        return rho


def base(residue, temperature, pH, value):

    if residue == "H":
        lamb = log10(10**(0.8-pH)/(10**(-get_pK_his(temperature))+10**(-pH))+10**(-0.1-get_pK_his(temperature))/(10**(-get_pK_his(temperature))+10**(-pH)))
        rho = log10(10**(0.83-pH)/(10**(-get_pK_his(temperature))+10**(-pH))+10**(0.14-get_pK_his(temperature))/(10**(-get_pK_his(temperature))+10**(-pH)))
    elif residue == "D":
        lamb = log10(10**(0.69-pH)/(10**(-get_pK_asp(temperature))+10**(-pH))+10**(0.1-get_pK_asp(temperature))/(10**(-get_pK_asp(temperature))+10**(-pH)))
        rho = log10(10**(0.6-pH)/(10**(-get_pK_asp(temperature))+10**(-pH))+10**(-0.18-get_pK_asp(temperature))/(10**(-get_pK_asp(temperature))+10**(-pH)))
    elif residue == "E":
        lamb = log10(10**(0.24-pH)/(10**(-get_pK_glu(temperature))+10**(-pH))+10**(-0.11-get_pK_glu(temperature))/(10**(-get_pK_glu(temperature))+10**(-pH)))
        rho = log10(10**(-0.39-pH)/(10**(-get_pK_glu(temperature))+10**(-pH))+10**(-0.15-get_pK_glu(temperature))/(10**(-get_pK_glu(temperature))+10**(-pH)))
    else:
        lamb = para[residue][2]
        rho = para[residue][3]

    if value == "lamb":
        return lamb
    elif value == "rho":
        return rho
