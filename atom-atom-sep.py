#!/usr.bin/env python
# Calculate the distribution of atom-atom separations
# relative to the sigma - defined diameter
# Input: atom dump files for cold film
# Calculate: mean and standard deviation of atom-atom separation
#            normalized to 2^(1/6)*sigma for each type pair
#            as a function of atom type.
# Output: List of normalized mean and stddev atom-atom separation for types 3 to 14
#
# R.B. Stephens Mar 2023
#
# Written for Python3
import os, sys, math
from numpy import *

# define measurement region
mz_min = 10.0
mz_max = 25.0
size_z = abs(mz_max-mz_min)

# atom-atom separation correction = sigma * 2 ^ (1/6)
d_cor  = 2**(1/6) # - correction for distance to potential minimum

box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

try:
  inatom = sys.argv[1]
except:
  print("Usage:", sys.argv[0], "(1) inatom")
  sys.exit(1)

if os.path.exists(inatom):
    mydirect = os.path.dirname(inatom)
    myfilename = os.path.basename(inatom)
else:
  print("Usage:", sys.argv[0], "(1) inatom")
  sys.exit(1)
mystring = myfilename.split('.')
if len(mystring)>0:
  basename = mystring[0]
  
# open and start reading atom file
readfile = mydirect + "/" + basename + ".film7"
outfile = "Reading atom file: " + readfile + "\n"
print(outfile)
iatom = open( readfile, 'r') # open file for reading

line = iatom.readline()
# check that first lines of atomfile has the correct word
if not line:
  print("missing first line in atomfile"); exit()
iwords = line.split()
if not iwords[0]=='LAMMPS':
  print("missing atomfile"); exit()

while 1:
  line = iatom.readline()
  if not line:
    print("atom and box parameters loop")
    print("missing line in atom and box param - atomfile"); exit()
  iwords = line.split()
  if len(iwords)==2:
    if iwords[1]=='atoms':
      at_no = int(iwords[0])
      print('total number of atoms in file  = %i' %(at_no))
      atm_x = zeros(at_no+1)
      atm_y = zeros(at_no+1)
      atm_z = zeros(at_no+1)
      atm_t = zeros(at_no+1,dtype=int)
      atm_n_cont = zeros(at_no+1,dtype=int)

  elif len(iwords)==3:
    if iwords[2]=='types':
      #     set up arrays for number of each type, and total and average contact number for each type
      at_type_no  = int(iwords[0]) # number of types
      atm_dia     = zeros((at_type_no+1,at_type_no+1)) # interatomic separation for each pair of types
      print('number of atom types = %i' %(at_type_no))
  elif len(iwords)>=4:
    if iwords[2]=='xlo':
    # get box x-dim
      box_x[0] = float(iwords[0])
      box_x[1] = float(iwords[1])
      size_x = abs(box_x[1]-box_x[0])
      mod_x = size_x/2.0
    elif iwords[2]=='ylo':
    # get box y-dim
      box_y[0] = float(iwords[0])
      box_y[1] = float(iwords[1])
      size_y = abs(box_y[1]-box_y[0])
      mod_y = size_y/2.0
    elif iwords[0]=='PairIJ':
      # start to read potentials
      while 1: # exits at beginning of Atoms section
      # read atom-atom sigma, put corrected value in array to get effective diameter
        line = iatom.readline()
        if not line:
          print("missing line reading sigmas in atomfile"); exit()
        iwords = line.split()
        if len(iwords)>0:
          if iwords[0]=='Atoms':
            print("finished sigmas - found atoms keyword")
            break
          atm0 = int(iwords[0])
          atm1 = int(iwords[1])
          atm_dia[atm0,atm1] = float(iwords[3])*d_cor
          atm_dia[atm1,atm0] = float(iwords[3])*d_cor
      break # go to next loop to read atom positions

# beginning of read atom position loop
iatm = -1
solid_vol = 0
while 1: #
  line = iatom.readline()
  if not line:
    print("missing line in atom positions - atomfile"); exit()
  iwords = line.split()
  if len(iwords)>0:
    if iwords[0]=='Velocities':
      break  # finished reading atom positions
    # read atom type, position
    # (I can't use given atom label - some numbers might have been skipped)
    iwords = line.split()
    if len(iwords)>5:
      iatm += 1
      atm_z[iatm] = float(iwords[4])
      atm_x[iatm] = float(iwords[2])
      atm_y[iatm] = float(iwords[3])
      atm_t[iatm] = int(iwords[1])

# total number of nearby pairs in the volume
ccount = 30*at_no
# calculate the separations and contacts for each nearby atom pair
norm_dist = zeros((at_no+1,ccount))
touch_norm_dist = zeros((at_no+1,ccount))

# list the contact number for each atom - the list sorted by type
# determine the normalized separation from all atoms within 1.5 of atom dia.
# create list and count of separations binned by atom type
# count number of separations <=1.0
# take mean and standard dev of list,

cont0 = zeros(at_no+1,dtype=int) # no of neighbors (separations less than 1.5)
cont1 = zeros(at_no+1,dtype=int) # no of contacts (separations less than 1.0)
atm_pr_cnt = zeros((at_no+1,at_no+1),dtype=int) # no of contacts of each type
natm_t = zeros(at_no+1,dtype=int) # no of atoms of each type
atm_no = iatm
# step through all the atoms in the characterized region
for a1ndx in range(atm_no+1):
  if ((atm_z[a1ndx]>= mz_min) and (atm_z[a1ndx]<= mz_max)):
     natm_t[atm_t[a1ndx]] += 1
     # look through all the other atoms to find the first's neighbors
     for a2ndx in range(atm_no+1):
       if a2ndx == a1ndx: continue
       sep_x = abs(atm_x[a1ndx]-atm_x[a2ndx])
       if (sep_x > mod_x):
         sep_x = sep_x - mod_x
       sep_y = abs(atm_y[a1ndx]-atm_y[a2ndx])
       if (sep_y > mod_y):
         sep_y = sep_y - mod_y
       sep_z = abs(atm_z[a1ndx]-atm_z[a2ndx])
       dist = sqrt(sep_x**2 + sep_y**2 + sep_z**2)
       act_sep = atm_dia[atm_t[a1ndx],atm_t[a2ndx]]
       testdist = dist/act_sep
       # count all the seps < 1.5 and put the values into a list  - both separated by type
       if testdist <=1.5:
         cont0[int(atm_t[a1ndx])] += 1
         norm_dist[atm_t[a1ndx],cont0[atm_t[a1ndx]]] = testdist
         # count all the seps <=1.0 - call those contacts - classify by atom2 type
         if testdist <= 1.0:
           cont1[atm_t[a1ndx]] += 1
           touch_norm_dist[atm_t[a1ndx],cont1[atm_t[a1ndx]]] = testdist
           atm_pr_cnt[atm_t[a1ndx],atm_t[a2ndx]] += 1
          
# calculate statistics for atom-atom separation for each atom type
atm_avg_sep = zeros(at_type_no+1)
atm_std_sep = zeros(at_type_no+1)
atm_avg_touch = zeros(at_type_no+1)
atm_avg_touch_sep = zeros(at_type_no+1)
atm_avg_touch_std = zeros(at_type_no+1)
atm_avg_neigh = zeros(at_type_no+1)
for andx in range (3,at_type_no+1):
  atm_avg_touch[andx] = cont1[andx]/natm_t[andx]
  atm_avg_neigh[andx] = cont0[andx]/natm_t[andx]
  atm_avg_sep[andx] = average(norm_dist[andx,0:cont0[andx]]) # average separation for atom of this type with all other atoms
  atm_std_sep[andx] = std(norm_dist[andx,0:cont0[andx]]) # std dev of separations for atom of this type with all other atoms
  atm_avg_touch_sep[andx] = average(touch_norm_dist[andx,0:cont1[andx]])
  atm_avg_touch_std[andx] = std(touch_norm_dist[andx,0:cont1[andx]])
  

# create contacts file, output results
fileword = myfilename.split('.')
outfile = "st2-" + fileword[0] + ".txt"
print(outfile)
out_stat = open(outfile, 'w')
# file header
out_stat.write('Data analysis for file: \n')
out_stat.write(basename)
out_stat.write(' \n')
out_stat.write('Atom-atom separation stats in selected region z= %g to %g \n' %(mz_min,mz_max))
out_stat.write('The separation and std-dev normalized to minimum in the atom-atom separation potential for each pair \n')
out_stat.write('Atom type, type count, separation, std dev, contact no, contact avg, std dev, neighbor no \n ')

for tndx in range(3,at_type_no+1):
  out_stat.write(' %7i, %7i, %10.5g, %10.5g, %10.5g, %10.5g, %10.5g, %10.5g \n' %(tndx,natm_t[tndx], atm_avg_sep[tndx], atm_std_sep[tndx], atm_avg_touch[tndx], atm_avg_touch_sep[tndx], atm_avg_touch_std[tndx], atm_avg_neigh[tndx]))
  
# matrix of atom-atom contact types
out_stat.write('Atom contact type matrix \n ')
outfile = " -- , "
for tndx in range(3,at_type_no+1):
  outfile = outfile + str(tndx) + ", "
outfile = outfile + " \n"
out_stat.write(outfile)
for a1ndx in range(3,at_type_no+1):
  outfile = str(a1ndx) + ", "
  for a2ndx in range(3,at_type_no+1):
    outfile = outfile + str(atm_pr_cnt[a1ndx,a2ndx]) + ", "
  outfile = outfile + " \n"
  out_stat.write(outfile)

print("Finished Calculations")
