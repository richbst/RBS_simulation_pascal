#!/usr.bin/env python

# Calculate density within height limits from an atom file -
# Input:       atomfile describing the set of atoms in this simulation
# Calculation: Average density vs distance from growing surface
# Output:      print calculated density
#
# R.B. Stephens October 2023
#
# Written for Python3
import sys, math, os
from numpy import *

try:
  inatom = sys.argv[1]
except:
  print("Usage:", sys.argv[0], "(1) inatom ")
  sys.exit(1)

# open file for reading the simulation parameters
iatom = open( inatom, 'r')

# create tag for output
infile = os.path.basename(inatom)
outtag = infile.split('.')
print (outtag[0])

box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

vol_coef = 3.1415926535/6.0
vol_cor  = sqrt(2) # correcting for distance to potential minimum

# read atom positions line by line
# the indices of each parameter are:
i_no = 0
i_typ = 1
i_x = 2
i_y = 3
i_z = 4
# the number of atom types will be read from file
ntype = 0

# calculate density within these bounds
avg_min = 10
avg_max = 25

# read first lines of atomfile to determine simulation parameters
#      (it skips atom masses - those are all set to 1)
#      (only read potential parameters between like sizes)
line = iatom.readline()
if not line:
  print("missing line in atomfile"); exit()
iwords = line.split()
if not iwords[0]=='LAMMPS':
  print("missing atomfile"); exit()
readpar = 1
while readpar:
  line = iatom.readline()
  iwords = line.split()
  if len(iwords)>2:
    if iwords[2]=='types':
      ntype = int(iwords[0])+1
      print('Number of atom types = %i' % (ntype))
      atm_dia  = zeros(ntype+1)
      atm_vol  = zeros(ntype+1)
      atm_no   = zeros(ntype+1)
    elif iwords[2]=='xlo':
      box_x[0] = float(iwords[0])
      box_x[1] = float(iwords[1])
      print('x_dim = (%.2f,%.2f)' % (box_x[0],box_x[1]))
    elif iwords[2]=='ylo':
      box_y[0] = float(iwords[0])
      box_y[1] = float(iwords[1])
      print('y_dim = (%.2f,%.2f)' % (box_y[0],box_y[1]))
    elif iwords[2]=='zlo':
      box_z[0] = float(iwords[0])
      box_z[1] = float(iwords[1])
      print('z_dim = (%.2f,%.2f)' % (box_z[0],box_z[1]))
    elif iwords[0]=='PairIJ':
      atm_typ_cnt = 1
      while 1:
        line = iatom.readline()
        if not line: print ('Error in PairIJ'); exit()
        iwords = line.split()
        if (len(iwords)>0):
          if (iwords[0]=='Atoms'): readpar = 0; print ('Found atom list');  break
          if (int(iwords[0])==atm_typ_cnt)&(int(iwords[1])==atm_typ_cnt):
            atm_dia[atm_typ_cnt] = float(iwords[3])
            atm_vol[atm_typ_cnt] = vol_coef*atm_dia[atm_typ_cnt]*atm_dia[atm_typ_cnt]*atm_dia[atm_typ_cnt]*vol_cor
            atm_typ_cnt += 1
            
atm_tot = 0
readatm = 1
while readatm:
  line = iatom.readline()
  if not line: readatm = 0; break
  iwords = line.split()
  if len(iwords)>=5:
    if(float(iwords[i_z])>=avg_min)&(float(iwords[i_z])<=avg_max):
      atm_no[int(iwords[i_typ])] += 1
      atm_tot += 1
  
print ('Finished atom list')
print ('Atom types = %i' % (atm_typ_cnt))
print ('Atom Total in region = %i' %(atm_tot))
for i in range (0,atm_typ_cnt,1):
  print ('Type(%i) = %i' % (i, atm_no[i]))
meas_vol = (avg_max-avg_min)*(box_x[1]-box_x[0])*(box_y[1]-box_y[0])
print ('Measured volume = %.2f' % (meas_vol))
tot_vol = 0
for i in range (0,ntype+1,1):
  tot_vol += atm_vol[i]*atm_no[i]
avg_dens = tot_vol/meas_vol
print (outtag)
print ('volume density between %i and %i is %.4g' %(avg_min, avg_max, avg_dens))
                    
                    
                                                
