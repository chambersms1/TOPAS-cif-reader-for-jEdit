# -*- coding: utf-8 -*-

___version___ = 1.1
"""
Created on 2023/05/05
Update: added extra tabs to the output file (2023/06/15)

@author: Matthew Chambers

Now produces a text file with correct string formatting
Ultimately want to make it compatible with XInsert in jedit
Some minor issues that can be adjusted by user but need to be sorted out:
    For ID, the sign of charge is in the wrong place for TOPAS
Perhaps make the script make "0.3333" into "=1/3;: 0"
Possibly make another version that accounts for U ADPs
Possibly update name so that it follows the "common name" e.g. "rutile". Not all .cif files use common name
"""

#import os
#import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("filename",help = "Input .cif file to be read")

args = parser.parse_args()
with open(args.filename, 'r') as readfile:
# 'cif' is all of the lines in the file that is read
# Each i in cif[i] corresponds to a line in the file
    cif = readfile.readlines()

# Define 'Atom' and 'Cell' classes as they have multiple properties
# Case sensitive!
class AtomCif(object):
    def __init__(self):
        self.name = ''
        self.type = ''
        self.id = 0
        self.OS = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occ = 0.0
        self.B = 0.0    # Isotopic atomic displacement parameter

class Cell(object):
    def __init__(self):
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.V = 0
        self.alpha = 0
        self.beta = 0
        self.gamma = 0
        self.formula = ''
        self.space_group = ''

# Script does an exhaustive search for lines that have certain strings at the beginning
# When it finds a line that has the correct terms, it splits the line so that the desired parameter from that line is extracted
atoms = []
for i in range(0, len(cif)):
    cif_line = cif[i]
    L0 = cif_line.split("'")
    L1 = cif_line.split()
    if str(cif_line[:28]) == "_chemical_formula_structural": #Some cifs do not use "'" for the formula, so needed to add this convoluted mess 
        if len(L0) == 1:
            L2 = L0[0].split("_chemical_formula_structural") #This splits it into list with \n in L2[1], so added extra parameter L3
            L3 = L2[1].replace("\n","")
            Cell.formula = L3.replace(" ","") #Removes the spaces between atoms
        else:
            Cell.formula = L0[1].replace(" ","")    # This is the "odd one out"; only line that needs to be split with "'"
    if str(cif_line[:14]) == "_cell_length_a":
        Cell.a = L1[1]
    if str(cif_line[:14]) == "_cell_length_b":
        Cell.b = L1[1]
    if str(cif_line[:14]) == "_cell_length_c":
       Cell.c = L1[1]
    if str(cif_line[:17]) == "_cell_angle_alpha":
        Cell.alpha = L1[1]
    if str(cif_line[:16]) == "_cell_angle_beta":
        Cell.beta = L1[1]
    if str(cif_line[:17]) == "_cell_angle_gamma":
        Cell.gamma = L1[1]
    if str(cif_line[:12]) == "_cell_volume":
        Cell.V = L1[1]
    if str(cif_line[:25]) == "_space_group_name_H-M_alt":
        Cell.space_group = L0[1].replace(" ","") # This also needs to have the lines split differently
    if str(cif_line[:30]) == "_symmetry_space_group_name_H-M":
        Cell.space_group = L0[1].replace(" ","") # This is the old format for the space group label; the executable in the TOPAS folder works with this format
    if str(cif_line[:20]) =="_atom_site_occupancy":
        for j in range((i+1), len(cif)):
            cif_line = cif[j]
            L = cif_line.split()
            if len(L) == 9:
                atom = AtomCif()
                atom.name = L[0]
                atom.id = L[1]
                atom.type = re.split('(\d+)',atom.id)[0] #Use this to get the element without the OS. Splits it from the first number, making a list of ['element', 'x', 'sign']
                atom.OS = re.split('(\d+)',atom.id)[2]+re.split('(\d+)',atom.id)[1] #Use this to rearrange the '+'/'-' with the number in OS; + or - will be the third entry in the list and the number of OS is the second entry
                atom.x = L[4]
                atom.y = L[5]
                atom.z = L[6]
                atom.B = L[7]
                atom.occ = L[8]
                atoms.append(atom)
                
with open("tempcif.txt",'w') as cif_temp:
            
#Added .split("(")[0] to the variables to remove the parentheses in the errors        
    cif_temp.write("\tstr")
    cif_temp.write('\n\t\tphase_name "'+Cell.formula+'"')
    cif_temp.write("\n\t\ta "+Cell.a.split("(")[0])
    cif_temp.write("\n\t\tb "+Cell.b.split("(")[0])
    cif_temp.write("\n\t\tc "+Cell.c.split("(")[0])
    cif_temp.write("\n\t\tal "+Cell.alpha.split("(")[0])
    cif_temp.write("\n\t\tbe "+Cell.beta.split("(")[0])
    cif_temp.write("\n\t\tga "+Cell.gamma.split("(")[0])
    cif_temp.write("\n\t\tvolume "+Cell.V.split("(")[0])
    cif_temp.write('\n\t\tspace_group "'+Cell.space_group+'"')
    #Use for loop to print the different atom qualities of each atom in one line
    for i in range(len(atoms)):
        cif_temp.write('{:10}'.format("\n\t\tsite "+atoms[i].name)+'{:10}'.format("\tx "+atoms[i].x.split("(")[0])+'{:10}'.format("\ty "+atoms[i].y.split("(")[0])+'{:10}'.format("\tz "+atoms[i].z.split("(")[0])+'{:10}'.format("\tocc "+atoms[i].type+atoms[i].OS)+'{:10}'.format(atoms[i].occ.split("(")[0])+'{:10}'.format("\tbeq "+atoms[i].B.split("(")[0]))

#Output for atomic sites without string formatting
#cif_temp.write("\nsite "+atoms[i].name+"\tx "+atoms[i].x.split("(")[0]+"\ty "+atoms[i].y.split("(")[0]+"\tz "+atoms[i].z.split("(")[0]+"\tocc "+atoms[i].id+" "+atoms[i].occ.split("(")[0]+"\tbeq "+atoms[i].B.split("(")[0])