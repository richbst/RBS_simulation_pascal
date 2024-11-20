#!/usr.bin/env python
# Calculate the distribution of relative sphere sizes of a grown film
# both the number at a distance and the number within a distance
# normalized by the number of atoms studied
#
# R.B. Stephens Oct 2024
# Written for Python3

import sys, math
from numpy import *

try:
  inatom = sys.argv[1];
except:
  print("Usage:", sys.argv[0], "(1) inatom")
  sys.exit(1)

print("start")

# box size - assume pbc for x and y
box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

# calculate distribution for all atoms in z = 15 to z = 20 for distance up to z = 8
z_min = 15.0
z_max = 20.0
r_max = 8.0
r_bin = 0.1
bin_max = int(r_max/r_bin + 0.5)

# constants used in atomfile
vol_coef = 3.1415926535/6.0
vol_cor  = sqrt(2) # correcting for distance to potential minimum

iatom = open( inatom, 'r') # open file for reading
print ("file opened to read")
file_words = inatom.split('/')
name_length = len(file_words)
print_name = file_words[name_length-1]
print (print_name)

# read first lines of atomfile
# determine number of types and sizes of each
line = iatom.readline()
if not line:
  print("missing line in atomfile"); exit()
iwords = line.split()
if not iwords[0]=='LAMMPS':
  print("missing atomfile"); exit()
readatm = 1
while readatm:
  line = iatom.readline()
  iwords = line.split()
  if len(iwords)==2:
    if iwords[1]=='atoms':
      # get atom number - set up atom spec array
      at_no = int(iwords[0])
      print ('got atom number = %i' % (at_no))
      atm_spec = zeros((at_no+1,6))
  if len(iwords)>2:
    if iwords[2]=='types':
      # get simulation dimensions
      print ("got space dimensions")
      at_type_no = int(iwords[0])
      atm_dia  = zeros(at_type_no+1)
      atm_vol  = zeros(at_type_no+1)
    elif iwords[2]=='xlo':
      box_x[0] = float(iwords[0])
      box_x[1] = float(iwords[1])
      x_size = (box_x[1]-box_x[0])
    elif iwords[2]=='ylo':
      box_y[0] = float(iwords[0])
      box_y[1] = float(iwords[1])
      y_size = (box_y[1]-box_y[0])
    elif iwords[0]=='Atoms':
      atm_cnt = 0
      atm_mid = 0
      print ("starting to read atom positions")
      while 1:
        line = iatom.readline()
        if not line:
          print ("file end")
          break
        iwords = line.split()
        if len(iwords)>0:
          if iwords[0]=='Velocities':
            readatm = 0
            print ("finished atom positions")
            break
          zz_pos = float(iwords[4])
          zz_type = int(iwords[1])
          if (zz_type>=3):
            atm_cnt+=1
            if (zz_pos>z_min) and (zz_pos<z_max):
              atm_mid+= 1
            atm_spec[atm_cnt,1] = atm_mid     # within calc region
            atm_spec[atm_cnt,2] = zz_type     # type
            atm_spec[atm_cnt,3] = iwords[2]   # x coord
            atm_spec[atm_cnt,4] = iwords[3]   # y coord
            atm_spec[atm_cnt,5] = zz_pos      # z coord
iatom.close();
print ("atom data file closed")
atm_tot = atm_cnt # number of atoms counting only types 3 to 14
mid_tot = atm_mid # number of atoms within calculation region
atm_sep = zeros((mid_tot+2,atm_tot+2,3))
out_string = ( 'data read: total atoms = %i, total within z = %i' % (atm_tot, mid_tot))
print (out_string)
# find all relevant separations
for jatm in range(1, atm_tot, 1):
  jmid = int(atm_spec[jatm,1])
  if jmid==0 : continue
  xj = atm_spec[jatm,3]
  yj = atm_spec[jatm,4]
  zj = atm_spec[jatm,5]
  for iatm in range(0, atm_tot, 1):
    if jatm==iatm : continue
    xi = atm_spec[iatm,3]
    yi = atm_spec[iatm,4]
    zi = atm_spec[iatm,5]
    del_x = abs(xi - xj)
    del_y = abs(yi - yj)
    if del_x>(x_size/2.0):
      del_x = x_size - del_x
    if del_y>(y_size/2.0):
      del_y = y_size - del_y
    del_z = abs(zi - zj)
    del_r = sqrt(del_x*del_x + del_y*del_y + del_z*del_z)
    if del_r<=r_max:
      atm_sep[jmid,iatm,0] = jatm
      atm_sep[jmid,iatm,1] = del_r
      typ_dif = int(abs(atm_spec[iatm,2] - atm_spec[jatm,2]))
      atm_sep[jmid,iatm,2] = typ_dif
    else:
      atm_sep[jmid,iatm,1] = -1
      atm_sep[jmid,iatm,2] = -1
print ("separations calculated")
#
# For number of atoms at given separation
# bin separations by type no diff
type_dif_max = 11
sep_distrib_s = zeros((type_dif_max+1,bin_max+1))
sep_scale = zeros(bin_max+1)
for indx in range(0,bin_max+1,1):
  sep_scale[indx] = indx*r_bin
for jatm in range(1, atm_mid, 1):
  for iatm in range(1, atm_tot, 1):
    separation = atm_sep[jatm,iatm,1]
    if separation==0: continue
    sep_indx = int(separation/r_bin +0.5)
    typ_dif = int(atm_sep[jatm,iatm,2])
    sep_distrib_s[typ_dif,sep_indx]+=1 # adds a count for this type diff and separation
print ("distribution calculated")
# correct count for number of atoms examined and  r^2 dependent sampling volume
rnorm_distrib_s = zeros((type_dif_max+1, bin_max+1))
countmax_s = 0
for rndx in range(1, bin_max, 1):
  corr = rndx*rndx*mid_tot
  count = 2*float(sep_distrib_s[0,rndx]/corr)
  rnorm_distrib_s[0,rndx] = count
  if count>countmax_s: countmax_s = count
  for dndx in range(1, type_dif_max, 1):
    count = float(sep_distrib_s[dndx,rndx]/corr)
    rnorm_distrib_s[dndx,rndx] = count
    if count>countmax_s: countmax_s = count
y_axis_max_s = countmax_s
print ("distance binning calculated")
#
# For number of atoms within given separation
# calculate distribution of all atoms with same type diff within bndx
sep_distrib_v = zeros((type_dif_max+1,bin_max+1))
for bndx in range(1,bin_max,1):
  for indx in range(1,bndx,1):
    for tndx in range(1,type_dif_max,1):
      sep_distrib_v[tndx,bndx]+=sep_distrib_s[tndx,indx]
# correct count for number of atoms examined and r^3 dependent sampling volume
rnorm_distrib_v = zeros((type_dif_max+1, bin_max+1))
countmax_v = 0
for rndx in range(1, bin_max, 1):
  corr = rndx*rndx*mid_tot
  count = 2*float(sep_distrib_v[0,rndx]/corr)
  rnorm_distrib_v[0,rndx] = count
  if count>countmax_v: countmax_v = count
  for dndx in range(0, type_dif_max, 1):
    count = float(sep_distrib_v[dndx,rndx]/corr)
    rnorm_distrib_v[dndx,rndx] = count
    if count>countmax_v: countmax_v = count
y_axis_max_v = countmax_v
print ("volume binning calculated")
# create text file of distribution
out_dist = open("diff_dist.txt", 'w')
# create header
out_string = ('Radial distribution distribution of atom size differences. \n')
out_dist.write(out_string)
out_string = ('From file:, ')
out_string = out_string + print_name + (', total atoms =, %i, mid atoms =, %i \n' % (atm_tot, mid_tot))
out_dist.write(out_string)
# write x-axis
out_string = (' , ')
for xndx in range (1, bin_max, 1):
  out_string = out_string + ('%5.2g, ' % (sep_scale[xndx]))
out_string = out_string + (' \n')
out_dist.write(out_string)
# write data number at separation
for typ_dif in range(0, type_dif_max, 1):
  out_string = ('%i, ' % (typ_dif))
  for bindx in range(1, bin_max, 1):
    out_string = out_string + ('%i, ' % (sep_distrib_s[typ_dif,bindx]))
  out_string = out_string + ('\n')
  out_dist.write(out_string)
out_string = ('\n\n\n')
print (out_string)
# write data number within separation
for typ_dif in range(0, type_dif_max, 1):
  out_string = ('%i, ' % (typ_dif))
  for bindx in range(1, bin_max, 1):
    out_string = out_string + ('%i, ' % (sep_distrib_v[typ_dif,bindx]))
  out_string = out_string + ('\n')
  out_dist.write(out_string)
out_dist.close()
# analyze curve height
y = zeros(bin_max)
# print graph of distribution at distance for each type diff number
import matplotlib.pyplot as plt
from matplotlib import cm
x = sep_scale
for indx in range(0,type_dif_max,1):
  y = rnorm_distrib_s[indx,:]
  label_str = ('diff %i' % (indx))
  plt.plot(x,y,label=label_str)
# label graph and axes
text_kwargs = dict(ha='left', va='center', fontsize=8, color='black', font='Helvetica', backgroundcolor='white')
x_pos = -0.12*r_max
y_pos = -0.1*y_axis_max_s
plt.text(x_pos,y_pos,print_name,**text_kwargs)
plt.xlabel('Binned at separation (units = %3.2f sigma)' % (r_bin))
plt.ylabel('Atom number')
plt.axis([0, r_max, 0, y_axis_max_s])
plt.title('Atom distribution at separation for size diff')
plt.figlegend()
plt.savefig('distrib_s.pdf')
plt.show()
# print graph of distribution within distance for each type diff number
import matplotlib.pyplot as plt
from matplotlib import cm
x = sep_scale
for indx in range(0,type_dif_max,1):
  y = rnorm_distrib_v[indx,:]
  label_str = ('diff %i' % (indx))
  plt.plot(x,y,label=label_str)
# label graph and axes
text_kwargs = dict(ha='left', va='center', fontsize=8, color='black', font='Helvetica', backgroundcolor='white')
x_pos = -0.12*r_max
y_pos = -0.1*y_axis_max_v
plt.text(x_pos,y_pos,print_name,**text_kwargs)
plt.xlabel('Binned within separation (units = %3.2f sigma)' % (r_bin))
plt.ylabel('Atom number')
plt.axis([0, r_max, 0, y_axis_max_v])
plt.title('Atom distribution within separation for size diff')
plt.figlegend()
plt.savefig('distrib_v.pdf')
plt.show()
print ("finished")
