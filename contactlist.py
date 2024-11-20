#!/usr.bin/env python
# Calculate the average number of contacts for each atom size
# This is meant to use atom dump files for films with no kinetic energy
#
# Inputs: atomfile describing the set of atoms in this simulation
#
# Calculate: ratio separation between each pair of atoms
#               to position of potential minimum for that pair
#            for each atom type sum the number for which ratio is <= 1
#
# Output: contact_no.txt. list of average number of contacts for each atom type, and number of that type
#            - the list containing atoms within the standard 10-25 height
#            - the contacts including atoms outside those boundaries
#  ????      - should I get number of contacts between pairs of types?
#  ????      - distribution of contact lengths - and normalized as fraction of equilibrium length
#
# R.B. Stephens Feb 2023
#
# Written for Python3
import os, sys, math
from numpy import *

# define measurement region
mz_min = 10.0
mz_max = 25.0

# 'contact overreach (fractional distance beyond equilibrium to claim contact)
atm_reach = 1.001

# atom-atom separation correction = sigma * 2 ^ (1/6)
d_cor  = 2**(1/6) # - correction for distance to potential minimum


try:
  inatom = sys.argv[1]
except:
  print("Usage:", sys.argv[0], "(1) inatom")
  sys.exit(1)
  
myfilename ="xxx"
if os.path.exists(inatom):
    myfilename = os.path.basename(inatom)

iatom = open( inatom, 'r') # open file for reading

box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

# read first lines of atomfile
# determine number of types and sizes of each
line = iatom.readline()
if not line:
  print("missing line in atomfile"); exit()
iwords = line.split()
if not iwords[0]=='LAMMPS':
  print("missing atomfile"); exit()

while 1:
  line = iatom.readline()
  if not line:
    print("missing line in atomfile"); exit()
  iwords = line.split()
  if len(iwords)>1: # there is a blank line at the beginning and end of each section
    if iwords[1]=='atoms':
      # use the number of atoms to set up atom info arrays
      #     for position, type, separation,
      #     atom-atom contact existence, and contact number
      at_no = int(iwords[0])
      atm_x  = zeros(at_no+1)
      atm_y  = zeros(at_no+1)
      atm_z  = zeros(at_no+1)
      atm_t  = zeros(at_no+1, dtype = int)
      print('total number of atoms  = %i' %(at_no))
    
    elif len(iwords)==3:
      if iwords[2]=='types':
        #     set up arrays for number of each type, and total and average contact number for each type
        at_type_no  = int(iwords[0]) # number of types
        atm_dia     = zeros((at_type_no+1,at_type_no+1)) # interatomic separation for each type pair
        t_atom_no    = zeros(at_type_no+1, dtype = int) # number of each type
        c_atm_t_no  = zeros(at_type_no+1) # number of contacts for each type
        t_con_avg   = zeros(at_type_no+1, dtype = float)  # average number of contacts for each type
        c_type_list = zeros((at_type_no+1,20))
        print('number of atom types = %i' %(at_type_no))
    
    elif len(iwords)==4:
      if iwords[2]=='xlo':
      # get box x-dim
        box_x[0] = float(iwords[0])
        box_x[1] = float(iwords[1])
        size_x = abs(box_x[1]-box_x[0])/2.0
        print("read box x values")
      elif iwords[2]=='ylo':
      # get box y-dim
        box_y[0] = float(iwords[0])
        box_y[1] = float(iwords[1])
        size_y = abs(box_y[1]-box_y[0])/2.0
        print("read box y values")
      elif iwords[0]=='PairIJ':
        break
        
# start read potentials loop
while 1: # exits at beginning of Atoms section
# read atom-atom sigma, put corrected value in array to get effective diameter
  line = iatom.readline()
  if not line:
    print("missing line in atomfile"); exit()
  iwords = line.split()
  if len(iwords)>0:
    if iwords[0]=='Atoms':
      break  # go to next loop to read atom positions
    atm0 = int(iwords[0])
    atm1 = int(iwords[1])
    atm_dia[atm0,atm1] = float(iwords[3])*d_cor
    atm_dia[atm1,atm0] = float(iwords[3])*d_cor
  
# beginning of read atom position loop
iatm = 0
while 1: #
  line = iatom.readline()
  if not line:
    print("missing line in atomfile"); exit()
  iwords = line.split()
  if len(iwords)>0:
    if iwords[0]=='Velocities':
      break  # finished reading atom positions
    # read atom type, position
    # (I can't use given atom label - some numbers might have been skipped)
    iwords = line.split()
    if len(iwords)>0:
      iatm += 1
      atm_x[iatm] = float(iwords[2])
      atm_y[iatm] = float(iwords[3])
      atm_z[iatm] = float(iwords[4])
      atm_t[iatm] = int(iwords[1])
    
atom_no = iatm
# atom_no might be smaller than at_no if some atoms went missing? Anyway trust the direct count
c_no = 0
c_t_indx = zeros(at_type_no+1, dtype = int) # number of contacts binned by type
c_t_length = zeros((at_type_no+1, 4000))
c_t_fractn = zeros((at_type_no+1, 4000))
# file each contact length under a list for that atom type - assume less than 20 contacts per atom, 100 atoms of each type.
for iatn1 in range(0, atom_no, 1):
  if ((atm_z[iatn1]>=mz_min) & (atm_z[iatn1]<=mz_max)):
  # only measure atoms inside the measurement region
    seltype = atm_t[iatn1]
    # print(seltype)
    t_atom_no[seltype] += 1
  # count the number of each type
    for iatn2 in range(1,atom_no,1):
    # find all atoms that are 'contacting' this atom - whether or not in the measurement region
    # use modulo operator for horizontal directions because of periodic b.c.
      xxx = abs(atm_x[iatn1]-atm_x[iatn2])
      yyy = abs(atm_y[iatn1]-atm_y[iatn2])
      zzz = abs(atm_z[iatn1]-atm_z[iatn2])
      if xxx > size_x:
        xxx -= size_x
      if yyy > size_y:
        yyy -= size_y
      # print(xxx, yyy, zzz)
      atm_sep = sqrt(xxx**2 + yyy**2 + zzz**2)
      atm1_type = atm_t[iatn1]
      atm2_type = atm_t[iatn2]
      atmatm_sep = atm_dia[atm1_type,atm2_type]
      atm_sep_fr = atm_sep/atmatm_sep
      # print(atm_sep, atmatm_sep)
      if atm_sep_fr<= atm_reach:
      # count them if 'contacting'
        c_t_indx[atm1_type] += 1 # counts total contact number for each type
        c_t_length[atm1_type,c_t_indx[atm1_type]] = atm_sep
        c_t_fractn[atm1_type,c_t_indx[atm1_type]] = atm_sep_fr

# calculate average contact number for each type
for ityp in range(3, at_type_no+1, 1):
  t_con_avg[ityp] = c_t_indx[ityp]/t_atom_no[ityp]

# calculate distribution - mean and standard deviation - sorted by type
# from numpy
cont_length = zeros(at_type_no+1)
cont_stdev  = zeros(at_type_no+1)
for ityp in range(1,at_type_no+1,1):
  cont_length[ityp] = mean(c_t_length[ityp,1:c_t_indx[ityp]]) # have to restrict to actual entries - don't accept zeros
  cont_stdev[ityp]  = std(c_t_length[ityp,1:c_t_indx[ityp]])/cont_length[ityp]

# output results - create contacts file
fileword = myfilename.split('.')
outfile = "contacts.txt"
if len(fileword)>0:
   outfile = "c-" + fileword[0] + ".txt"
out_contacts = open(outfile, 'w')
# file header
out_contacts.write('Data analysis for file: \n')
out_contacts.write(myfilename)
out_contacts.write(' \n')
out_contacts.write('Average contact number per atom overall and for each type \n')
out_contacts.write('Type, Atom No., Avg contact no., length, stdev/length \n ')
for indx in range(3, at_type_no+1, 1):
  out_contacts.write(' %5i, %5i, %7.3g, %7.3g, %7.3g \n' %(indx, t_atom_no[indx], t_con_avg[indx], cont_length[indx], cont_stdev[indx]))

out_contacts.close()
print("Finished Calculations")
