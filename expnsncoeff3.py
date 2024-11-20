#!/usr.bin/env python
# Calculate the film density at a set of temperatures
# Inputs: atom dump files for films as they are cooled.
#         ending in .film<num> with <num> from 1 to 6
#         (only put the first file as input - program figures out where to find the rest
# Calculate: density in selected region for each file
#
# Output: Temp / density list
#
# R.B. Stephens Feb 2023
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

temp_dens = zeros((8,3)) #temp, dens, box width

# temps used in LAMMPS script
temp_dens[0,0] = 0.21
temp_dens[1,0] = 0.15
temp_dens[2,0] = 0.1
temp_dens[3,0] = 0.03
temp_dens[4,0] = 0.01
temp_dens[5,0] = 0.003
temp_dens[6,0] = 0.001
temp_dens[7,0] = 0.0003

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
  
# loop through 8 atom files
fndx = -1
while fndx < 7:
  fndx += 1
  readfile = mydirect + "/" + basename + ".film" + str(fndx)
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
        print('total number of atoms in file[%i]  = %i' %(fndx, at_no))
    elif len(iwords)==3:
      if iwords[2]=='types':
        #     set up arrays for number of each type, and total and average contact number for each type
        at_type_no  = int(iwords[0]) # number of types
        atm_dia     = zeros(at_type_no+1) # interatomic separation for each type paired with itself
        a_vol  = zeros(at_no+1) # volume of each type
        print('number of atom types = %i' %(at_type_no))
    elif len(iwords)>=4:
      if iwords[2]=='xlo':
      # get box x-dim
        box_x[0] = float(iwords[0])
        box_x[1] = float(iwords[1])
        size_x = abs(box_x[1]-box_x[0])
        temp_dens[fndx,2] = size_x
        size_x = size_x/2.0
        print("read box x values")
      elif iwords[2]=='ylo':
      # get box y-dim
        box_y[0] = float(iwords[0])
        box_y[1] = float(iwords[1])
        size_y = abs(box_y[1]-box_y[0])/2.0
        print("read box y values")
      elif iwords[0]=='PairIJ':
        # calculate volume of measured region
        reg_vol = size_x*size_y*size_z*4
        # start read potentials loop
        while 1: # exits at beginning of Atoms section
        # read atom-atom sigma, put corrected value in array to get effective diameter
          line = iatom.readline()
          if not line:
            print("missing line reading potentials in atomfile"); exit()
          iwords = line.split()
          if len(iwords)>0:
            if iwords[0]=='Atoms':
              print("found Atoms keyword")
              break
            atm0 = int(iwords[0])
            atm1 = int(iwords[1])
            if atm0==atm1:
                atm_dia[atm0] = float(iwords[3])*d_cor
        for vndx in range (3, at_type_no+1, 1):
            a_vol[vndx] = pi/6*atm_dia[vndx]**3
        break # go to next loop to read atom positions

  # beginning of read atom position loop
  iatm = 0
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
        atm_z = float(iwords[4])
        atm_t = int(iwords[1])
        if ((atm_z>=mz_min) & (atm_z<=mz_max)):
        # only count atoms inside the measurement region
           solid_vol += a_vol[atm_t]
           # add the volume of each atom to solid volume

  temp_dens[fndx,1] = solid_vol/reg_vol
  # then go to next atom position file

# output results - create contacts file
fileword = myfilename.split('.')
outfile = "d3-" + fileword[0] + ".txt"
print(outfile)
out_density = open(outfile, 'w')
# file header
out_density.write('Data analysis for file: \n')
out_density.write(basename)
out_density.write(' \n')
out_density.write('Density in selected region \n')
out_density.write('Temp, Volume density, box width \n ')

for fndx in range(0, 8, 1):
  out_density.write(' %10.5g, %10.5g, %10.5g \n' %(temp_dens[fndx,0], temp_dens[fndx,1], temp_dens[fndx,2]))

out_density.close()

print("Finished Calculations")
